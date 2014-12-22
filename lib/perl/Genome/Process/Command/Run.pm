package Genome::Process::Command::Run;

use strict;
use warnings FATAL => 'all';
use Genome;
use Try::Tiny qw(try catch);
use Data::Dump qw(pp);
use Genome::Utility::Email;
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

    my %inputs = %$inputs;
    while (my ($name, $value) = each %inputs) {
        if (ref($value) eq 'HASH') {
            $inputs->{$name} = convert_hash_to_obj($value);
        } elsif (ref($value) eq 'ARRAY' &&
                 scalar(@{$value}) &&
                 ref($value->[0]) eq 'HASH') {
            $inputs->{$name} = [map {convert_hash_to_obj($_)} @{$value}];
        }
    }

    return $inputs
}

sub convert_hash_to_obj {
    my $hash = shift;

    my $class = $hash->{'class'};
    my $obj = $class->get($hash->{'id'});
    unless (defined($obj)) {
        die sprintf("Couldn't convert hash to class: %s", pp($hash));
    }
    return $obj;
}

sub execute {
    my $self = shift;

    $self->process->write_environment_file;

    $self->status_message("Reading in workflow from file: %s",
        $self->process->workflow_file);
    my $dag = Genome::WorkflowBuilder::DAG->from_xml_filename(
        $self->process->workflow_file);

    $dag->log_dir($self->process->log_directory);
    $dag->name($self->process->workflow_name);

    $self->update_status('Running');
    my $inputs = $self->get_workflow_inputs;
    $self->status_message("Executing workflow with inputs: %s", pp($inputs));

    my $error;
    my $outputs = try {
        $dag->execute(%$inputs);
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

sub schedule {
    my $self = shift;

    my $job_id = $self->_bsub_in_pend_state;
    $self->_ensure_bkill_if_rollback($job_id);
    $self->_ensure_bresume_if_commit($job_id);
}

sub _bsub_in_pend_state {
    my $self = shift;

    my $cmd = ['annotate-log', 'genome', 'process', 'run',
        '--process', $self->process->id];
    my @bsub_args = (
        cmd => $cmd,
        email => Genome::Utility::Email::construct_address(),
        err_file => $self->workflow_process_error_log,
        log_file => $self->workflow_process_out_log,
        hold_job => 1,
        project => $self->process->workflow_name,
        queue => $ENV{GENOME_LSF_QUEUE_BUILD_WORKFLOW},
        send_job_report => 1,
        never_rerunnable => 1,
    );

    my $job_id = try {
        Genome::Sys::LSF::bsub::run(['bsub'], @bsub_args);
    } catch {
        die "Failed to launch bsub:\n$_\n";
    };
    $self->process->update_status('Scheduled');
    return $job_id;
}

sub _ensure_bkill_if_rollback {
    my $self = shift;
    my $job_id = shift;

    my $bsub_undo = sub {
        $self->status_message("Killing LSF job (%s) for process (%s).",
            $job_id, $self->process->id);
        system("bkill $job_id");
    };

    my $lsf_change = UR::Context::Transaction->log_change($self, 'UR::Value', $job_id, 'external_change', $bsub_undo);
    unless ($lsf_change) {
        die $self->error_message("Failed to record LSF job submission ($job_id).");
    }
}

sub _ensure_bresume_if_commit {
    my $self = shift;
    my $job_id = shift;

    my $commit_observer = UR::Context->process->add_observer(
        aspect => 'commit',
        callback => sub {
            my $bresume_output = `bresume $job_id`; chomp $bresume_output;
            $self->error_message($bresume_output) unless ( $bresume_output =~ /^Job <$job_id> is being resumed$/ );
        },
    );
    unless ($commit_observer) {
        $self->error_message("Failed to add commit observer to resume LSF job ($job_id).");
    }
}


1;
