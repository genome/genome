package Genome::Model::Build::Command::Start;

use strict;
use warnings;

use Genome;

use Genome::Utility::Instrumentation;
use Time::HiRes;

use Data::Dumper 'Dumper';
use Regexp::Common;

use Try::Tiny qw(try catch);

class Genome::Model::Build::Command::Start {
    is => 'Command::V2',
    roles => 'Genome::Role::CommandWithColor',
    doc => "Create and start a build.",
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            doc => 'Model(s) to build. Resolved from command line via text string.',
            shell_args_position => 1,
        },
    ],
    has_optional => [
        data_directory => { },
        force => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Force a new build even if existing builds are running.',
        },
        skip_succeeded => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Prevent launching a new build if the model has a status of "Succeeded".',
        },
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
            is_output => 1,
        },
        max_builds => {
            is => 'Integer',
        },
        unstartable_ok => {
            is => 'Boolean',
            default => 0,
        },
        _builds_started => {
            is => 'Integer',
            default => 0,
        },
        _create_params => {
            is => 'Hash',
            default => {},
        },
        _start_params => {
            is => 'Hash',
            default => {},
        },
    ],

};

sub sub_command_sort_position { 1 }

sub help_synopsis {
    return <<EOS;
genome model build start 1234

genome model build start somename

EOS
}

sub help_detail {
    return <<EOS;
Make a new build for the specified model, and initiate execution of the build processes.

Builds with a defined workflow will run asynchronously.  Simple builds will run immediately
and this command will wait for them to finish.
EOS
}

sub execute {
    my $self = shift;

    $self->_create_params->{data_directory} = $self->data_directory if ($self->data_directory);

    my @models = $self->models;
    for my $model (@models) {
        if ($self->max_builds && $self->_builds_started >= $self->max_builds){
            $self->status_message("Already started max builds $self->_builds_started, quitting");
            last; 
        }
        $self->_total_command_count($self->_total_command_count + 1);
        if(!$self->force) {
            my @existing_builds = $model->builds(status => ['Running', 'Scheduled']);
            if (@existing_builds) {
                $self->append_error($model->__display_name__, "Model already has running or scheduled builds. Use the '--force' option to override this and start a new build.");
                next;
            }
        }
        if ($self->skip_succeeded and $model->status eq 'Succeeded') {
            my $msg = 'Skipping succeeded model ' . $model->__display_name__;
            $self->status_message($self->_color($msg, 'cyan'));
            next;
        }

        Genome::Utility::Instrumentation::timer('command.model.build.start', sub {
            $self->create_and_start_build($model);
        });
    }

    $self->display_builds_started();
    $self->display_command_summary_report();


    return !scalar(keys %{$self->_command_errors});
}

sub create_and_start_build {
    my $self = shift;
    my $model = shift;

    $self->status_message("Trying to start #" . ($self->_builds_started + 1) . ': ' . $model->__display_name__ . "...");

    my $outer_transaction = UR::Context::Transaction->begin();
    $model->build_requested(0);
    my $anp = $model->analysis_project;
    my $guard;
    $guard = $anp->set_env if $anp;

    my $create_transaction = UR::Context::Transaction->begin();
    my $build = try {
        my $build = Genome::Model::Build->create(model_id => $model->id, %{$self->_create_params})
            or die 'failed to create build';
        $self->add_build($build);
        $create_transaction->commit()
            or die 'failed to commit';
        return $build;
    }
    catch {
        $create_transaction->rollback();
        $outer_transaction->rollback();
        $self->append_error($model->__display_name__, $_);
        return;
    };
    unless ($build) {
        return;
    }

    my $start_transaction = UR::Context::Transaction->begin();
    my $build_started = try {
        my $rv = $build->start(%{$self->_start_params});
        $start_transaction->commit() or die "Cannot commit 'build start' transaction";
        return $rv;
    }
    catch {
        $start_transaction->rollback();
        $self->append_error($model->__display_name__, 'Build (' . $build->id . ') ' . $_);
        $build->status('Unstartable');
        $build->add_note(
            header_text => 'apipe_cron_status',
            body_text => $_,
            auto_truncate_body_text => 1,
        );
        return;
    };
    unless ($outer_transaction->commit) {
        die 'transaction failed to commit';
    }

    if (not $build_started) {
        return ($build->status eq 'Unstartable' && $self->unstartable_ok);
    }

    $self->_builds_started($self->_builds_started + 1);
    my $msg = "Successfully started build (" . $build->__display_name__ . ").";
    $self->status_message($self->_color($msg, 'green'));
    return 1;
}


sub display_builds_started {
    my $self = shift;
    my @builds = $self->builds;

    if (@builds) {
        $self->status_message("Build IDs: " . join(' ', map { $_->id } @builds));
    }
    return 1;
}


1;
