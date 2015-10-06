package Genome::Model::ReferenceVariation::Command::QualityControl;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceVariation::Command::QualityControl {
    is => 'Genome::Model::ReferenceVariation::Command::Base',
    doc => 'Perform quality control analysis on the alignments for the build.',
};

sub _result_accessor {
    return 'output_result';
}

sub _command_class {
    return 'Genome::Qc::Run';
}

sub _params_for_command {
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
