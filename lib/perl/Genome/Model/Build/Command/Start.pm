package Genome::Model::Build::Command::Start;

use strict;
use warnings;

use Genome;

use Genome::Utility::Instrumentation;
use Time::HiRes;

use Data::Dumper 'Dumper';
use Regexp::Common;

class Genome::Model::Build::Command::Start {
    is => 'Genome::Command::Base',
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
        job_dispatch => {
            doc => 'dispatch specification: an LSF queue or "inline"',
        },
        server_dispatch => {
            doc => 'dispatch specification: an LSF queue or "inline"',
        },
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
# default values for dispatching will be either -s workflow -j apipe
# or come from the processing profile if available as a param

genome model build start somename -s workflow -j apipe
# run the server in the workflow queue, and jobs in the apipe queue

genome model build start somename -s inline -j inline
# run the server inline, and the jobs inline

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
    $self->_start_params->{job_dispatch} = $self->job_dispatch if ($self->job_dispatch);
    $self->_start_params->{server_dispatch} = $self->server_dispatch if ($self->server_dispatch);

    my @models = $self->models;
    for my $model (@models) {
        if ($self->max_builds && $self->_builds_started >= $self->max_builds){
            $self->status_message("Already started max builds $self->_builds_started, quitting");
            last; 
        }
        $self->_total_command_count($self->_total_command_count + 1);
        if (!$self->force && ($model->builds_with_status('Running') or $model->builds_with_status('Scheduled'))) {
            $self->append_error($model->__display_name__, "Model already has running or scheduled builds. Use the '--force' option to override this and start a new build.");
            next;
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
    my $create_transaction = UR::Context::Transaction->begin();
    my $build = eval {
        my $build = Genome::Model::Build->create(model_id => $model->id, %{$self->_create_params});
        unless ($build) {
            die $self->error_message($model->__display_name__, "Failed to create new build.");
        }
        return $build;
    };

    if ($build and $create_transaction->commit) {
        # Record newly created build so other tools can access them.
        # TODO: should possibly be part of the object class
        $self->add_build($build);

        my $start_transaction = UR::Context::Transaction->begin();
        my $build_started = eval { $build->start(%{$self->_start_params}) };
        if ($start_transaction->commit) {
            if ($build_started) {
                $self->_builds_started($self->_builds_started + 1);
                my $msg = "Successfully started build (" . $build->__display_name__ . ").";
                $self->status_message($self->_color($msg, 'green'));
            }
            else {
                if ($build->status eq 'Unstartable') {
                    $self->append_error($model->__display_name__, 'Build (' . $build->id . ') created but Unstartable, review build\'s notes.');
                }
                elsif ($@) {
                    $self->append_error($model->__display_name__, 'Build (' . $build->id . ') ' . $@);
                }
                else {
                    $self->append_error($model->__display_name__, 'Build (' . $build->id . ') not started but unable to parse error, review console output.');
                }
            }
        }
        else {
            # If we couldn't commit after trying to start then something blocked us from even committing that the build was Unstartable.
            my $transaction_error = $start_transaction->error_message();
            $start_transaction->rollback;
            my $error_message = 'Failed to commit build start, rolling back to build creation.';
            if ($transaction_error) {
                $error_message .= " Transaction error: $transaction_error";
            }
            $self->append_error($model->__display_name__, $error_message);
            $build->status('Unstartable');
            $build->add_note(
                header_text => 'apipe_cron_status',
                body_text => $error_message,
                auto_truncate_body_text => 1,
            );
            $build->model->build_requested(0);
        }
    }
    else {
        if ($@) {
            print "Exception is $@";
            $self->append_error($model->__display_name__, $@);
        }
        else {
            $self->append_error($model->__display_name__, 'Build not created but unable to parse error, review console output.');
        }
        $create_transaction->rollback;
    }
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

