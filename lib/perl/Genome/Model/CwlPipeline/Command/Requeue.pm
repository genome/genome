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

    my $count = 0;
    BUILD: for my $build ($self->builds) {

        unless ($build->status eq 'Failed') {
            $self->error_message('Cannot restart build %s with status %s. Only Failed builds may be restarted.', $build->__display_name__, $build->status);
            next BUILD;
        }

        unless($build->run_by eq Genome::Config::get('builder_run_as_user')) {
            $self->error_message('Can only requeue builds for the configured builder user (%s). Restart your own builds directly with `genome model cwl-pipeline restart`.', Genome::Config::get('builder_run_as_user'));
            next BUILD;
        }

        $build->rebuild_requested(1);

        my @note_params = (
            header_text => 'Rebuild Requested',
        );
        if ($self->reason) {
            push @note_params,
                body_text => $self->reason;
        }
        $build->add_note(@note_params);

        $count++;
    }

    $self->status_message('Requeued %s failed build(s).', $count);

    return 1;
}

1;
