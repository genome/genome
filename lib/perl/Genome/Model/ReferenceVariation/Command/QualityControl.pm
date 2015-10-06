package Genome::Model::ReferenceVariation::Command::QualityControl;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceVariation::Command::QualityControl {
    is => 'Command::V2',
    has_input => [
        build => {
            is => 'Genome::Model::Build::ReferenceVariation',
            id_by => 'build_id',
            doc => 'build for which to run quality control',
            is_output => 1,
        },
    ],
    doc => 'Perform quality control analysis on the alignments for the build.',
};

sub shortcut {
    my $self = shift;

    my $cmd = $self->_quality_control_command;

    my $retval = $cmd->shortcut;
    if ($retval and $cmd->output_result) {
        $self->status_message('Found existing result: %s', $cmd->output_result->id);
    }

    return $retval;
}

sub execute {
    my $self = shift;

    my $cmd = $self->_quality_control_command;
    my $retval = $cmd->execute;

    if ($retval and $cmd->output_result) {
        $self->status_message('Generated result: %s', $cmd->output_result->id);
    } else {
        die $self->error_message('Failed to produce qc result.');
    }

    return $retval;
}

sub _quality_control_command {
    my $self = shift;

    my @params = $self->_params_for_quality_control;

    my $cmd = Genome::Qc::Run->create(@params);
    die $self->error_messaage('Failed to create qc command!') unless $cmd;

    return $cmd;
}

sub _params_for_quality_control {
    my $self = shift;
    my $build = $self->build;

    my $result_users = Genome::SoftwareResult::User->user_hash_for_build($build);
    $result_users->{qc_result} = $build;

    my @params = (
        alignment_result => $build->merged_alignment_result,
        config_name => $build->model->qc_config,
        result_users => $result_users,
    );

    return @params;
}

1;
