package Genome::Model::Build::Command::Purge;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::Purge {
    is => 'Command::V2',
    has         => [
        builds => {
            is                  => 'Genome::Model::Build',
            is_many             => 1,
            shell_args_position => 1,
            doc                 => "Builds whose data to purge",
            require_user_verify => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    my @builds = $self->builds;

    for my $build (@builds) {
        $build->purge or $self->fatal_message('Failed to purge build %s', $build->__display_name__);
    }
}

1;
