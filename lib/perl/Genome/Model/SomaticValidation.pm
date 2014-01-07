package Genome::Model::SomaticValidation;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation {
    is  => 'Genome::ModelDeprecated',
    has_param_optional => [
        alignment_strategy => {
            is => 'Text',
            is_many => 0,
            doc => 'Strategy to be used to align',
        },
        snv_detection_strategy => {
            is => "Text",
            is_many => 0,
            doc => "Strategy to be used to detect snvs.",
        },
        indel_detection_strategy => {
            is => "Text",
            is_many => 0,
            doc => "Strategy to be used to detect indels.",
        },
        sv_detection_strategy => {
            is => "Text",
            is_many => 0,
            doc => "Strategy to be used to detect svs.",
        },
        cnv_detection_strategy => {
            is => "Text",
            is_many => 0,
            doc => "Strategy to be used to detect cnvs.",
        },
        identify_dnp_proportion => {
            is => 'Number',
            doc => 'The proportion of reads supporting a DNP to make the call',
        },
        minimum_coverage => {
            is => 'Number', 
            doc => 'minimum coverage to call a site (in process-validation step)',
        },
        output_plot => {
            is => 'Boolean', 
            doc => 'include output plot in final results',
        },
        loh_version => {
            is => 'Text',
            doc => 'Version of LOH to use. (LOH detection is not performed if not specified.)',
        },
        loh_snv_detection_strategy => {
            is => 'Text',
            doc => 'The detection to run on the control aligned reads for determining LOH',
        },
        transcript_variant_annotator_version => {
            doc => 'Version of the "annotate transcript-variants" tool to run during the annotation step',
            is_optional => 1,
            default_value => Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version,
            valid_values => [ 0,1,2,3,4 ],
        },
        tiering_version => {
            is => 'Text',
            doc => 'version of tiering BED files to use (tiering is not performed if not specified)',
        },
        varscan_validation_version => {
            is => 'Text',
            doc => 'version of varscan to use in the post-variant detection validation processes',
        },

        #RefCov parameters
        refcov_wingspan_values => {
            is => 'Text',
            doc => 'Comma-delimited list of wingspans to use',
        },
        refcov_minimum_depths => {
            is => 'Text',
            doc => 'Comma-delimited list of depth levels to use',
        },
        refcov_minimum_base_quality => {
            is => 'Text',
            doc => 'Minimum base quality for consideration',
        },
        refcov_minimum_mapping_quality => {
            is => 'Text',
            doc => 'Minimum mapping quality for consideration',
        },
        refcov_merge_roi_regions => {
            is => 'Boolean',
            doc => 'Merge contiguous regions in the ROI set before analysis',
        },
        refcov_use_short_names => {
            is => 'Boolean',
            doc => 'Replace names in the BED file with short pre-generated names',
        },
        refcov_roi_track_name => {
            is => 'Text',
            valid_values => ['target_region','tiled_region'],
            doc => 'For multi-tracked BED files we define which track to use.',
        },
    ],
    has_input_optional_mutable => [
        snv_variant_list => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
        },
        indel_variant_list => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
        },
        sv_variant_list => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
        },
        reference_sequence_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
        },
        previously_discovered_variations_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'build of variants to screen out from consideration (such as from dbSNP)',
        },
        target_region_set => {
            is => 'Genome::FeatureList',
        },
        region_of_interest_set => {
            is => 'Genome::FeatureList',
        },
        design_set => {
            is => 'Genome::FeatureList',
        },
        tumor_sample => {
            is => 'Genome::Sample',
        },
        normal_sample => {
            is => 'Genome::Sample',
        },
        known_sites => {
            is => 'Genome::Model::Build::ImportedVariationList',
            is_many => 1,
            doc => 'Build[s] of known variants to use in when refining with GATK best practices.',
        },
    ],
    has_optional => [
        region_of_interest_set_name => {
            is => 'Text',
            via => 'region_of_interest_set',
            to => 'name',
        },
        target_region_set_name => {
            is => 'Text',
            via => 'target_region_set',
            to => 'name',
        },
        experimental_subject => {
            is => 'Genome::Sample',
            via => '__self__',
            to => 'tumor_sample'
        },
        control_subject => {
            is => 'Genome::Sample',
            via => '__self__',
            to => 'normal_sample'
        },

    ],
    has_transient_constant_optional => {
        sequencing_platform => {
            value => undef,
            doc => 'This can be removed once it has been removed from Genome::Model',
            is_deprecated => 1,
        },
    },
};

sub _validate_required_for_start_properties {
    my $self = shift;

    my @missing_required_properties;
    push @missing_required_properties, 'reference_sequence_build' unless ($self->reference_sequence_build);
    push @missing_required_properties, 'tumor_sample' unless ($self->tumor_sample);
    push @missing_required_properties, 'instrument_data' unless (scalar @{[ $self->instrument_data ]} );
    push @missing_required_properties, 'annotation_build' unless ($self->annotation_build);

    my $tag;
    if (@missing_required_properties) {
        $tag = UR::Object::Tag->create(
            type => 'error',
            properties => \@missing_required_properties,
            desc => 'missing required property',
        );
    }

    return $tag;
}

#limit compatible instrument data check to these samples
sub get_all_possible_samples {
    my $self = shift;

    my @result;
    push @result, $self->tumor_sample if $self->tumor_sample;
    push @result, $self->normal_sample if $self->normal_sample;

    return @result;
}

sub _resolve_workflow_for_build {
    my $self = shift;
    my $build = shift;

    my $operation = Workflow::Operation->create_from_xml(__FILE__ . '.xml');

    my $log_directory = $build->log_directory;
    $operation->log_dir($log_directory);
    $operation->name($build->workflow_name);

    return $operation;
}


sub map_workflow_inputs {
    my $self = shift;
    my $build = shift;

    my @inputs = ();

    # Verify the somatic model
    my $model = $build->model;
    unless ($model) {
        $self->error_message("Failed to get a model for this build!");
        die $self->error_message;
    }

    my $data_directory = $build->data_directory;
    unless ($data_directory) {
        $self->error_message("Failed to get a data_directory for this build!");
        die $self->error_message;
    }

    my $reference_sequence_build = $model->reference_sequence_build;
    unless($reference_sequence_build) {
        $self->error_message("Failed to get a reference sequence build for this model!");
        die $self->error_message;
    }
    my $reference_fasta = $reference_sequence_build->full_consensus_path('fa');
    unless(Genome::Sys->check_for_path_existence($reference_fasta)) {
        $self->error_message('Could not find reference FASTA for specified reference sequence.');
        die $self->error_message;
    }

    push @inputs,
        build_id => $build->id,
        transcript_variant_annotator_version => $build->processing_profile->transcript_variant_annotator_version,
        tumor_mode => 'tumor',
        normal_mode => 'normal',
        ;

    my %default_filenames = $self->default_filenames;
    for my $param (keys %default_filenames) {
        my $default_filename = $default_filenames{$param};
        push @inputs,
            $param => join('/', $data_directory, $default_filename);
    }

    push @inputs,
        minimum_coverage => (defined $self->minimum_coverage ? $self->minimum_coverage : 0),
        output_plot => (defined $self->output_plot ? $self->output_plot : 1),
        ;


    return @inputs;
}

sub default_filenames{
    my $self = shift;

    my %default_filenames = (
        targeted_snv_validation => 'validation/metrics/targeted.snvs.validation',
    );

    return %default_filenames;
}

sub default_profile {
    return Genome::ProcessingProfile::SomaticValidation->get(
        name => 'Feb 2013 default Somatic Validation Extension and Targeted Discovery');
}

sub default_single_bam_profile {
    return Genome::ProcessingProfile::SomaticValidation->get(
        name => 'Jun 2012 Single-Bam Validation (single-bam somatic)');
}

sub requires_pairing { return 1; }

1;
