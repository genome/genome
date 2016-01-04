package Genome::Process::Command::DiffBlessed;

use strict;
use warnings;
use Genome;

class Genome::Process::Command::DiffBlessed {
    is => 'Command::V2',
    roles => [qw/ Genome::Role::Comparable::Command::Diff Genome::Role::Comparable::Command::DiffBlessed /],
    has => [
        new_process => {
            is => 'Genome::Process',
            doc => "The process that you'd like to know if it has the correct output.",
            shell_args_position => 1,
        },
        process_name => {
            is => "String",
        },
    ],
};

sub help_brief : Overrides(Genome::Role::Comparable::Command::Diff, Genome::Role::Comparable::Command::DiffBlessed) {
    &Genome::Role::Comparable::Command::DiffBlessed::help_brief;
}

sub help_detail : Overrides(Genome::Role::Comparable::Command::Diff, Genome::Role::Comparable::Command::DiffBlessed) {
    &Genome::Role::Comparable::Command::DiffBlessed::help_detail;
}

sub new_object {
    my $self = shift;
    return $self->new_process;
}

sub blessed_object {
    my $self = shift;
    my $blessed_id = $self->blessed_id($self->process_name);
    return Genome::Process->get(id => $blessed_id);
}

1;

