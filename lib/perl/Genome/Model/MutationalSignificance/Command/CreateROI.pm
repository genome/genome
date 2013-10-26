package Genome::Model::MutationalSignificance::Command::CreateROI;

use strict;
use warnings;

use Genome;

class Genome::Model::MutationalSignificance::Command::CreateROI {
    is => ['Command::V2'],
    has_input => [
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation'
        },
        excluded_reference_sequence_patterns => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            doc => "Exclude transcripts on these reference sequences",
        },
        included_feature_type_patterns => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            doc => 'Include only entries that match one of these patterns',
        },
        include_flank => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Include the flanking regions of transcripts',
        },
        condense_feature_name => {
            is => 'Boolean',
            doc => 'Use only gene name as feature name',
            default_value => 1,
        },
        flank_size => {
            is => 'Integer',
            doc => 'Add this number of base pairs on each side of the feature', #to do: check this
            default_value => 0,
        },
        extra_rois => {
            is => 'Genome::FeatureList',
            is_many => 1,
            is_optional => 1,
        },
        include_ensembl_annot => {
            is => 'Boolean',
            default_value => 1,
        },
        regulome_bed => {
            is => 'Genome::FeatureList',
            is_optional => 1,
            doc => 'Bed file of regions with regulomedb scores',
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => '-R \'select[mem>32000] rusage[mem=32000]\' -M 32000000 ',
        },
    ],
    has_output => [
        roi_path => {
            is => 'String',
        },
    ],
};

sub execute {
    my $self = shift;
    my $new_feature_list;
    if ($self->include_ensembl_annot) {
        my @params;

        if ($self->excluded_reference_sequence_patterns) {
            push @params, excluded_reference_sequence_patterns => [$self->excluded_reference_sequence_patterns];
        }
        if ($self->included_feature_type_patterns) {
            push @params, included_feature_type_patterns => [$self->included_feature_type_patterns];
        }
        if ($self->condense_feature_name) {
            push @params, condense_feature_name => $self->condense_feature_name;
        }
        if ($self->flank_size && $self->flank_size > 0) {
            push @params, flank_size => $self->flank_size;
        }
        if ($self->include_flank) {
            push @params, include_flank => 1;
        }
        push @params, one_based => 1;

        my $feature_list = $self->annotation_build->get_or_create_roi_bed(@params);

        unless ($feature_list) {
            $self->error_message('ROI file not available from annotation build '.$self->annotation_build->id);
            return;
        }

        $self->status_message("Basing ROI on ".$feature_list->id);

        $self->roi_path($feature_list->file_path);

        my $new_name = $feature_list->name;
        my ($names, @files) = $self->collect_files_from_extra_rois($new_name);
        $new_name = $names;
        
        $new_feature_list = Genome::FeatureList->get(name => $new_name);

        unless ($new_feature_list) {
            $new_feature_list = $self->create_new_feature_list($new_name, $feature_list->subject, $feature_list->reference, $feature_list->format, @files);
            unless ($new_feature_list) {
                $self->error_message("Failed to create ROI file with extra ROIs");
                return;
            }
        }
        if ($self->regulome_bed) {
            my $filtered_name = join("_", $new_feature_list->name, "filtered_by_regulome_v1");
            my $filtered_list = Genome::FeatureList->get($filtered_name);
            unless ($filtered_list) {
                #transform to 0-based
                my $zero_based = $new_feature_list->processed_bed_file(
                    short_name => 0,
                );

                #filter
                my $filtered_out_zero_based = Genome::Sys->create_temp_file_path;
                my $rv = Genome::Model::Tools::RegulomeDb::ModifyRoisBasedOnScore->execute(
                    roi_list => $zero_based,
                    scored_regions => $self->regulome_bed->file_path,
                    output_file => $filtered_out_zero_based,
                    valid_scores => [qw(1 2)],
                );

                #convert back to 1-based
                my $filtered_out = Genome::FeatureList::transform_zero_to_one_based(
                    $filtered_out_zero_based,
                    $new_feature_list->is_multitracked,
                );

                my $file_content_hash = Genome::Sys->md5sum($filtered_out);
                my $filtered_feature_list = Genome::FeatureList->create(
                    name => $filtered_name,
                    format => $new_feature_list->format,
                    file_content_hash => $file_content_hash,
                    subject => $new_feature_list->subject,
                    reference => $new_feature_list->reference,
                    file_path => $filtered_out,
                    content_type => "roi",
                    description => "Feature list with extra rois filtered by regulome db",
                    source => "WUTGI",
                );

                $new_feature_list = $filtered_feature_list;
            }
        }
    }
    else {
        my ($names, @files) = $self->collect_files_from_extra_rois("rois");
        $new_feature_list = Genome::FeatureList->get(name => $names);
        unless($new_feature_list) {
            my @extras = $self->extra_rois;
            my $subject = $extras[0]->subject;
            my $reference = $extras[0]->reference;
            my $format = $extras[0]->format;
            $new_feature_list = $self->create_new_feature_list($names, $subject, $reference, $format, @files);
        }
    }
    
    $self->roi_path($new_feature_list->file_path);
    $self->status_message('Using ROI file: '.$self->roi_path);
    return 1;
}

sub create_new_feature_list {
    my $self = shift;
    my $name = shift;
    my $subject = shift;
    my $reference = shift;
    my $format = shift;
    my @files = @_;
    my $sorted_out = Genome::Sys->create_temp_file_path;
    my $rv = Genome::Model::Tools::Joinx::Sort->execute(
        input_files => [@files],
        unique => 1,
        output_file => $sorted_out,
    );
    my $file_content_hash = Genome::Sys->md5sum($sorted_out);


    my $feature_list = Genome::FeatureList->create(
        name => $name,
        format => $format,
        file_content_hash => $file_content_hash,
        subject => $subject,
        reference => $reference,
        file_path => $sorted_out,
        content_type => "roi",
        description => "Feature list with extra rois",
        source => "WUTGI",
    );
    return $feature_list;
}

sub collect_files_from_extra_rois {
    my $self = shift;
    my $name = shift;
    my @files;
    for my $extra_roi ($self->extra_rois) {
        my $roi_name = $extra_roi->name;
        $self->status_message("Adding roi $roi_name");
        $name .= "_$roi_name";
        push @files, $extra_roi->get_one_based_file;
    }
    return ($name, @files);
}

1;
