package Genome::Model::Build::Command::Diff;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::Diff {
    is => 'Genome::Interfaces::Comparable::Command::Diff',
    has => [
        blessed_build => {
            is => 'Genome::Model::Build',
            doc => "The build which is known to have the correct output.",
            shell_args_position => 2,
        },
        new_build => {
            is => "Genome::Model::Build",
            shell_args_position => 1,
        },
    ],
};

sub blessed_object {
    my $self = shift;
    return $self->blessed_build;
}

sub new_object {
    my $self = shift;
    return $self->new_build;
}
