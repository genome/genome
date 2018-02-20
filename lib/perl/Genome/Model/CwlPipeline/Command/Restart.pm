package Genome::Model::CwlPipeline::Command::Restart;

use strict;
use warnings;

use File::Spec;
use Genome;

class Genome::Model::CwlPipeline::Command::Restart {
    is => 'Command::V2',
    has => [
        builds => {
            is => 'Genome::Model::Build::CwlPipeline',
            is_many => 1,
            doc => 'The build(s) to restart',
        },
    ],
    doc => 'relaunch a failed build',
};

sub help_detail {
    return <<EOHELP
Restart a build process after a failure has occurred.  This will attempt to resume the in-progress workflow.
EOHELP
    ;
}

sub execute {
    my $self = shift;

    my @builds = $self->builds;
    my $build_count = scalar(@builds);
    my @errors;
    for my $build (@builds) {
        my $transaction = UR::Context::Transaction->begin();
        local $@;
        my $successful = eval { $self->_restart_build($build) };
        my $error = $@ // "";
        if ($successful and $transaction->commit) {
            $self->status_message(
                'Build (%s) restarted. An initialization email will be sent once the build begins running.',
                $build->__display_name__,
            );
        }
        else {
            $self->error_message(
                'Failed to restart build (%s): %s.',
                $build->__display_name__,
                $error
            );
            $transaction->rollback;
        }
    }

    return !scalar(@errors);
}

sub _restart_build {
    my $self = shift;
    my $build = shift;

    unless($build->status eq 'Failed') {
        die 'Can only restart failed builds.';
    }

    my $anp = $build->model->analysis_project;
    my $guard;
    $guard = $anp->set_env if $anp;

    $build->status('Scheduled');
    $build->date_completed(undef);

    my %params; #need to pass empty hash to launch
    my $xml = Genome::Sys->read_file(
        File::Spec->join($build->data_directory, 'build.xml')
    );

    return $build->_launch(\%params, $xml);
}

1;
