package Genome::Model::CwlPipeline::Command::Requeue;

use strict;
use warnings;

use File::Spec;
use Genome;

class Genome::Model::CwlPipeline::Command::Requeue {
    is => 'Command::V2',
    has => [
        builds => {
            is => 'Genome::Model::Build::CwlPipeline',
            is_many => 1,
            doc => 'The build(s) to queue for restart',
            shell_args_position => 1,
        },
        action => {
            is => 'Text',
            doc => 'what action to request for the build (or "0" to cancel an existing request")',
            valid_values => Genome::Model::CwlPipeline->__meta__->property(property_name => 'action_requested')->valid_values,
            default => 'restart',
        },
    ],
    has_optional => [
        reason => {
            is => 'Text',
            doc => 'a note on why the build(s) are being requeued, if desired',
        },
    ],
    doc => 'queue a reattempt for a failed build',
};

sub help_detail {
    return <<EOHELP
Queues a build to be restarted.  This is only useful for "production" builds.  For user-launched builds, use the `restart` command instead.
EOHELP
}

sub execute {
    my $self = shift;

    my $action = $self->action;

    my $count = 0;
    BUILD: for my $build ($self->builds) {

        if ($action eq 'restart') {
            if ($build->status ne 'Failed') {
                $self->error_message('Cannot restart build %s with status %s. Only Failed builds may be restarted.', $build->__display_name__, $build->status);
                next BUILD;
            }

            unless($build->run_by eq Genome::Config::get('builder_run_as_user')) {
                $self->error_message('Can only requeue builds for the configured builder user (%s). Restart your own builds directly with `genome model cwl-pipeline restart`.', Genome::Config::get('builder_run_as_user'));
                next BUILD;
            }
        } else {
            unless($build->run_by eq Genome::Config::get('builder_run_as_user')) {
                $self->error_message('Can only %s builds for the configured builder user (%s). For your own builds, use `genome model build %s` directly.', $action, Genome::Config::get('builder_run_as_user'), $action);
                next BUILD;
            }

            unless($build->model->analysis_project->created_by eq Genome::Sys->username or Genome::Sys->current_user_is_admin) {
                $self->error_message('Cannot %s builds for an analysis project owned by %s. (If this is your AnP, see `genome-analysis-project take-ownership`.', $action, $build->model->analysis_project->created_by);
                next BUILD;
            }
        }

        $build->action_requested($action);

        my @note_params = (
            header_text => 'Action Requested: ' . $action,
        );
        if ($self->reason) {
            push @note_params,
                body_text => $self->reason;
        }
        $build->add_note(@note_params);

        $count++;
    }

    $self->status_message('Requested %s for %s build(s).', $action, $count);

    return 1;
}

1;
