package Genome::Model::RnaSeq::Command::PicardRnaSeqMetrics;

use strict;
use warnings;

use Try::Tiny;

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
        picard_strand_specificity => {
            valid_values => Genome::Model::Tools::Picard::CollectRnaSeqMetrics->__meta__->property("strand_specificity")->valid_values,
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

sub should_skip {
    my $self = shift;

    my $build = $self->build;

    # Skip if no annotation build
    unless ( $build->annotation_build ) {
        $self->debug_message('Skipping PicardRnaSeqMetrics since annotation build is not defined');
        return 1;
    }

    # Skip if annotation build does not have required files
    my $rv;
    try {
        Genome::InstrumentData::AlignmentResult::Command::PicardRnaSeqMetrics->verify_annotation_build_has_required_files(
            $build->annotation_build, $build->reference_build,
        );
    }
    catch {
        $self->debug_message($_);
        $self->debug_message('Skipping PicardRnaSeqMetrics since annotation build is missing required files');
        $rv = 1;
    };

    return $rv;
}

sub shortcut {
    my $self = shift;

    return 1 if $self->should_skip;

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

    return 1 if $self->should_skip;

    my $build = $self->build;
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
        my $cmd = Genome::InstrumentData::AlignmentResult::Command::PicardRnaSeqMetrics->create(%params);
        unless($cmd and $cmd->execute) {
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

    my $result_users = Genome::SoftwareResult::User->user_hash_for_build($build);
    $result_users->{picard_rna_seq_metrics} = $build;

    my %params = (
        alignment_result_id => $alignment_result->id,
        picard_version => $self->picard_version,
        test_name => (Genome::Config::get('software_result_test_name') || undef),
        users => $result_users,
    );
    if ($self->picard_strand_specificity) {
        $params{picard_strand_specificity} = $self->picard_strand_specificity;
    } elsif ($build->model->picard_strand_specificity) {
        $params{picard_strand_specificity} = $build->model->picard_strand_specificity;
    }
    return %params;
}

sub link_result_to_build {
    my $self = shift;
    my $result = shift;

    my $build = $self->build;
    Genome::Sys->create_symlink($result->output_dir, $build->metrics_directory);

    $self->picard_rna_seq_metrics_result($result);

    return 1;
}


1;
