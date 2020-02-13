package Genome::Process::Command::SetStatus;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::Process::StatusEvent;

class Genome::Process::Command::SetStatus {
    is => ['Command::V2'],
    has => [
        process => {
            is => 'Genome::Process',
            shell_args_position => 1,
            doc => "The Genome::Process that you want to update the status on.",
        },
        status => {
            is => 'Text',
            shell_args_position => 2,
            valid_values => Genome::Process::StatusEvent->valid_status_values(),
            doc => 'The new status.',
        },
        exit_code => {
            is => 'Number',
            # We usually want to exit 1 so the step 'fails', but sometimes we want to
            # exit 0.
            default => '1',
            doc => 'The exit code for this command.',
        },
    ],
    doc => 'Sets the status of a Genome::Process.',
};

sub help_detail {
    my $self = shift;
    return <<EOP;
Used internally by the Genome::Process infrastructure to set the status of a process
while it is being run in a workflow.
EOP
}

sub shortcut {
    my $self = shift;
    return $self->execute();
}

sub execute {
    my $self = shift;

    $self->process->update_status($self->status);
    UR::Context->commit();
    exit($self->exit_code);
}


1;
