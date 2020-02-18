package Genome::Process::Command::Run;

use strict;
use warnings FATAL => 'all';
use Genome;
use Try::Tiny qw(try catch);
use Data::Dump qw(pp);
use Genome::Utility::Email;
use Genome::Utility::Inputs qw(decode);
use JSON qw(from_json);

class Genome::Process::Command::Run {
    is => ['Command::V2'],
    has => [
        process => {
            is => 'Genome::Process',
            shell_args_position => 1,
            doc => 'Process to run',
        },
        update_with_commit => {
            is => 'Boolean',
            default => 1,
            doc => 'Execute UR::Context->commit() after each status update',
        },
    ],
    doc => 'This command is not designed for users to run.',
};

sub help_synopsis {
    return <<EOS;
genome process run <process_id>
EOS
}

sub help_detail {
    return <<EOS;
This command is not designed for users to run.  It is the part of the
Genome::Process framework that allows processes to be run in a non-blocking
manner.
EOS
}

sub shortcut {
    return 0;
}

sub get_workflow_inputs {
    my $self = shift;

    my $json = Genome::Sys->read_file($self->process->inputs_file);
    my $inputs = from_json($json);

    return decode($inputs);
}

sub execute {
    my $self = shift;

    my $dag = $self->get_dag_for_execution;

    $self->update_status('Running');
    my $inputs = $self->get_workflow_inputs;
    $self->status_message("Executing workflow with inputs: %s", pp($inputs));

    my $error;
    my $outputs = try {
        $dag->execute(inputs => $inputs);
    } catch {
        $error = $_;
        $self->error_message("Workflow failed to complete: $error");
        $self->update_status('Crashed');
    };

    if ($error) {
        return 0;
    } else {
        $self->update_status('Succeeded');
        $self->status_message("Workflow completed with outputs: %s", pp($outputs));
        return $outputs;
    }
}

sub get_dag_for_execution {
    my $self = shift;

    $self->process->write_environment_file;

    $self->status_message("Reading in workflow from file: %s",
        $self->process->workflow_file);
    my $dag = Genome::WorkflowBuilder::DAG->from_xml_filename(
        $self->process->workflow_file);

    $dag->recursively_set_log_dir($self->process->log_directory);
    $dag->name($self->process->workflow_name);

    return $dag;
}

sub submit {
    my $self = shift;

    my $dag = $self->get_dag_for_execution;

    $self->update_status('Scheduled');

    my $inputs = $self->get_workflow_inputs;
    $self->status_message("Submitting workflow with inputs: %s", pp($inputs));

    my %env_copy = %ENV;

    my $commit_observer = Genome::Sys::CommitAction->create(
        on_commit => sub {
            local %ENV = %env_copy;
            my $wf_proxy = $dag->submit(inputs => $inputs, process => $self->process);
            $self->status_message("Successfully launched process (%s) and ".
                "submitted workflow (%s)", $self->process->id, $wf_proxy->url);
        },
    );
    unless ($commit_observer) {
        $self->error_message(sprintf "Failed to add commit observer to "
            ."submit workflow for process (%s).", $self->process->id);
    }
}

sub update_status {
    my ($self, $status) = @_;

    $self->process->update_status($status);
    if ($self->update_with_commit) {
        UR::Context->commit;
    }
}

sub workflow_process_error_log {
    my $self = shift;
    return File::Spec->join($self->process->log_directory, 'workflow.err');
}

sub workflow_process_out_log {
    my $self = shift;
    return File::Spec->join($self->process->log_directory, 'workflow.out');
}

1;
