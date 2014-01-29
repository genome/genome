package Genome::Model::RnaSeq::Command::Coverage;

use strict;
use warnings;

use Genome;
use version;

my $DEFAULT_LSF_RESOURCE = "-R 'select[mem>=32000] rusage[mem=32000]' -M 32000000";
my @DEFAULT_ANNOTATION_FILE_BASENAMES = qw/annotation rRNA rRNA_protein MT pseudogene/;
my $DEFAULT_MERGE_ANNOTATION_FEATURES = 'both';

class Genome::Model::RnaSeq::Command::Coverage {
    is => ['Command::V2'],
    has_input_output => [
        build_id => {},
    ],
    has => [
        build => { is => 'Genome::Model::Build::RnaSeq', id_by => 'build_id', },
        model => { via => 'build', },
    ],
    has_param => [
        lsf_resource => { default_value => $DEFAULT_LSF_RESOURCE },
    ],
    has_optional_output => [
        transcriptome_coverage_result_id => {
            is => 'Number',
            doc => 'ID of the result from running this step',
        },
    ],
    has_optional => [
        transcriptome_coverage_result => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged::TranscriptomeCoverage',
            id_by => 'transcriptome_coverage_result_id',
            doc => 'The result from running this step',
        },
    ],
};

sub sub_command_category { 'pipeline steps' }

sub sub_command_sort_position { 7 }

sub shortcut {
    my $self = shift;
    my $build = $self->build;
    
    my $pp = $build->processing_profile;
    if ($pp->transcriptome_coverage_annotation_file_basenames eq 'none') {
        $self->debug_message('The annotation file defined in the processing profile is \'none\'.  Transcriptome coverage will be skipped.');
        return 1;
    }
    my $annotation_build = $build->annotation_build;
    unless ($annotation_build) {
        $self->debug_message('Skipping TranscriptomeCoverage since annotation_build is not defined');
        return 1;
    }
    my $alignment_result = $build->alignment_result;
    if ($alignment_result->isa('Genome::InstrumentData::AlignmentResult::Merged')) {
        my %params = $self->params_for_result;
        my $result = Genome::InstrumentData::AlignmentResult::Merged::TranscriptomeCoverage->get_with_lock(%params);

        if($result) {
            $self->debug_message('Using existing result ' . $result->__display_name__);
            return $self->link_result_to_build($result);
        }
    }
    return;
}

sub execute {
    my $self = shift;

    my $build = $self->build;
    my $pp = $build->processing_profile;
    if ($pp->transcriptome_coverage_annotation_file_basenames eq 'none') {
        $self->debug_message('The annotation file defined in the processing profile is \'none\'.  Transcriptome coverage will be skipped.');
        return 1;
    }
    my $annotation_build = $build->annotation_build;
    unless ($annotation_build) {
        $self->debug_message('Skipping TranscriptomeCoverage since annotation_build is not defined');
        return 1;
    }
    my $alignment_result = $self->build->alignment_result;

    # Tophat v1.1.0 and later produces BAM output
    unless  (version->parse($alignment_result->aligner_version) >= version->parse('1.1.0')) {
        die('Coverage requires a BAM file produced by TopHat v1.1.0 or greater');
    }
    
    if ($alignment_result->isa('Genome::InstrumentData::AlignmentResult::Merged')) {
        my %params = (
            $self->params_for_result,
            log_directory => $build->log_directory,
        );
        my $result = Genome::InstrumentData::AlignmentResult::Merged::TranscriptomeCoverage->get_or_create(%params);
        $self->link_result_to_build($result);
    } else {
        # Run the old way...
        my %params =  (
            $self->params_for_result,
            coverage_directory => $build->coverage_directory,
            annotation_build => $build->annotation_build,
            reference_build => $build->reference_sequence_build,
        );
        #This is not a software result, yet...
        delete($params{test_name});
        unless (Genome::InstrumentData::AlignmentResult::Command::TranscriptomeCoverage->execute(%params)) {
            return;
        }
    }
    return 1;
}

sub default_annotation_file_basenames {
    my $class = shift;
    return \@DEFAULT_ANNOTATION_FILE_BASENAMES;
}

sub default_merge_annotation_features {
    my $class = shift;
    return $DEFAULT_MERGE_ANNOTATION_FEATURES;
}

sub params_for_result {
    my $self = shift;
    my $build = $self->build;
    my $pp = $build->processing_profile;

    my $alignment_result = $build->alignment_result;

    unless ($alignment_result) {
        die $self->error_message('No alignment result found for build: '. $build->id);
    }

    # Resolve MergeAnnotationFeatures
    my $merge_annotation_features;
    unless ($pp->transcript_coverage_merge_annotation_features) {
        $merge_annotation_features = $self->default_merge_annotation_features;
    } else {
        $merge_annotation_features = $pp->transcript_coverage_merge_annotation_features;
    }

    # Resolve AnnotationFileBasenames
    my $annotation_file_basenames;
    unless ($pp->transcriptome_coverage_annotation_file_basenames) {
        $annotation_file_basenames = $self->default_annotation_file_basenames;
    } else {
        my $string = $pp->transcriptome_coverage_annotation_file_basenames;
        my @annotation_file_basenames = split(',',$string);
        $annotation_file_basenames = \@annotation_file_basenames;
    }
    unless ($merge_annotation_features && $annotation_file_basenames) {
        die('Must define a default merge_annotation_features and/or annotation_file_basenames parameter!');
    }
    my $mask_reference_transcripts = undef;
    if ($pp->transcriptome_coverage_mask_reference_transcripts) {
        $mask_reference_transcripts = $pp->transcriptome_coverage_mask_reference_transcripts;
    }
    return (
        alignment_result_id => $alignment_result->id,
        annotation_file_basenames => $annotation_file_basenames,
        merge_annotation_features => $merge_annotation_features,
        mask_reference_transcripts => $mask_reference_transcripts,
        test_name => ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef),
    );
}

sub link_result_to_build {
    my $self = shift;
    my $result = shift;

    my $build = $self->build;
    my $label = join('_', 'transcriptome_coverage');
    Genome::Sys->create_symlink($result->output_dir, $build->coverage_directory);
    $result->add_user(label => $label, user => $build);

    $self->transcriptome_coverage_result($result);

    return 1;
}


1;

