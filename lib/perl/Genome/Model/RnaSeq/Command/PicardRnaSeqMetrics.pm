package Genome::Model::RnaSeq::Command::PicardRnaSeqMetrics;

use strict;
use warnings;

use Genome;

my $DEFAULT_LSF_RESOURCE = "-R 'select[mem>=8000] rusage[mem=8000]' -M 8000000";

class Genome::Model::RnaSeq::Command::PicardRnaSeqMetrics {
    is => ['Command::V2'],
    has_input_output => [
        build_id => {},
    ],
    has => [
        picard_version => {
            is_optional => 1,
            is_input => 1,
        },
        build => { is => 'Genome::Model::Build', id_by => 'build_id', },
    ],
    has_param => [
        lsf_resource => { default_value => $DEFAULT_LSF_RESOURCE },
    ],
    has_optional_output => [
        picard_rna_seq_metrics_result_id => {
            is => 'Number',
            doc => 'ID of the result from running this step',
        },
    ],
    has_optional => [
        picard_rna_seq_metrics_result => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged::PicardRnaSeqMetrics',
            id_by => 'picard_rna_seq_metrics_result_id',
            doc => 'The result from running this step',
        },
    ],
    doc => 'generate metrics on pipeline results'
};

sub sub_command_category { 'pipeline steps' }

sub sub_command_sort_position { 9 }

sub shortcut {
    my $self = shift;

    my $build = $self->build;
    my $alignment_result = $build->alignment_result;
    if ($alignment_result->isa('Genome::InstrumentData::AlignmentResult::Merged')) {
        my %params = $self->params_for_result;
        my $result = Genome::InstrumentData::AlignmentResult::Merged::PicardRnaSeqMetrics->get_with_lock(%params);
        if ($result) {
            $self->debug_message('Using existing result ' . $result->__display_name__);
            return $self->link_result_to_build($result);
        }
    }
    return;
}

sub execute {
    my $self = shift;

    my $build = $self->build;
    my $annotation_build = $build->annotation_build;
    unless ($annotation_build) {
        $self->debug_message('Skipping PicardRnaSeqMetrics since annotation_build is not defined');
        return 1;
    }
    my $alignment_result = $build->alignment_result;
    if ($alignment_result->isa('Genome::InstrumentData::AlignmentResult::Merged')) {
        my %params = (
            $self->params_for_result,
            log_directory => $build->log_directory,
        );
        my $result = Genome::InstrumentData::AlignmentResult::Merged::PicardRnaSeqMetrics->get_or_create(%params);
        $self->link_result_to_build($result);
    } else {
        # Run the old way...
        my %params = (
            $self->params_for_result,
            metrics_directory => $build->metrics_directory,
            annotation_build => $build->annotation_build,
            reference_build => $build->reference_sequence_build,
        );
        #This is not a software result, yet...
        delete($params{test_name});
        unless (Genome::InstrumentData::AlignmentResult::Command::PicardRnaSeqMetrics->execute(%params)) {
            return;
        }
    }
    return 1;
}

sub params_for_result {
    my $self = shift;
    
    my $build = $self->build;
    unless ($self->picard_version) {
        $self->picard_version($build->model->picard_version);
    }
    my $alignment_result = $build->alignment_result;

    unless ($alignment_result) {
        die $self->error_message('No alignment result found for build: '. $build->id);
    }

    return (
        alignment_result_id => $alignment_result->id,
        picard_version => $self->picard_version,
        test_name => ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef),
    );
}

sub link_result_to_build {
    my $self = shift;
    my $result = shift;

    my $build = $self->build;
    my $label = join('_', 'picard_rna_seq_metrics');
    Genome::Sys->create_symlink($result->output_dir, $build->metrics_directory);
    $result->add_user(label => $label, user => $build);

    $self->picard_rna_seq_metrics_result($result);

    return 1;
}


1;
