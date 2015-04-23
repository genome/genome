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
            valid_values => $Genome::Process::StatusEvent::VALID_STATUS_VALUES,
            doc => 'The new status.',
        }
    ],
    doc => 'Sets the status of a Genome::Process.',
};

sub help_detail {
    my $self = shift;
    return <<EOP;
Used internally by the Genome::Process infrastructure to set the status of a process
while it is being run in a workflow using PTero.
EOP
}

sub shortcut {
    my $self = shift;
    return $self->execute();
}

sub execute {
    my $self = shift;

    $self->process->update_status($self->status);
}


1;
