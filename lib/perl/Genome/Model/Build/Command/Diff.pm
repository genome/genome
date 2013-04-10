package Genome::Model::Build::Command::Diff;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::Diff {
    is => 'Genome::Model::Build::Command::Diff::Base',
    has => [
        blessed_build => {
            is => 'Genome::Model::Build',
            doc => "The build which is known to have the correct output.",
            shell_args_position => 2,
        },
    ],
};

sub get_blessed_build {
    my $self = shift;
    return $self->blessed_build;
}
