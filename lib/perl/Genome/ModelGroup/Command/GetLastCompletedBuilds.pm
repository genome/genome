package Genome::ModelGroup::Command::GetLastCompletedBuilds;

use strict;
use warnings;
use Genome;

class Genome::ModelGroup::Command::GetLastCompletedBuilds {
    is => 'Genome::Command::Base',
    doc => "Obtain all of the last_completed_builds from a model group. Not intended for command line use.",
    has_input => [
        model_group => {
            shell_args_position => 1,
            is => 'Genome::ModelGroup',
            doc => 'the model group from which to get builds', 
        },
    ],
    has_optional_output => [
        builds => {
            is_many => 1,
            is => 'Genome::Model::Build',
            doc => 'The builds that were resolved',
        },
    ],
};

sub execute {
    my $self = shift;
    my @builds;

    for my $model ($self->model_group->models) {
        my $build = $model->last_complete_build;
        if ($build) {
            push @builds, $build;
        }
        else {
            die $self->error_message("Found no last_completed_build for model " . $model->__display_name__);
        }
    }
    $self->builds(\@builds);

    $self->status_message(scalar(@builds) . " build(s) found in model group " . $self->model_group->__display_name__);

    return @builds;
};

1;

