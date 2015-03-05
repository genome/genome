package Genome::Process::Command::Diff;

use strict;
use warnings;

use Genome;

class Genome::Process::Command::Diff {
    is => 'Genome::Interfaces::Comparable::Command::Diff',
    has => [
        new_process => {
            is => 'Genome::Process',
            doc => "The process that you'd like to know if it has the correct output.",
            shell_args_position => 1,
        },
        blessed_process => {
            is => 'Genome::Process',
            shell_args_position => 2,
        },
    ],
};

sub blessed_object {
    my $self = shift;
    return $self->blessed_process;
}

sub new_object {
    my $self = shift;
    return $self->new_process;
}
