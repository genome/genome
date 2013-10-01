package Genome::Model::Build::Command::DelayArchiving;

use strict;
use warnings;
use Genome;

class Genome::Model::Build::Command::DelayArchiving {
    is => 'Command::V2',
    has => [
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
            shell_args_position => 1,
            doc => 'builds are marked so they cannot be archived',
        },
    ],
    has_optional => [
        Genome::Disk::Command::Allocation::DelayArchiving->_duration_property_for_commands()
    ],
};

sub help_detail {
    return 'Archiving of all allocations that are used by this build are delayed for 3 months' .
        'This includes allocations the build itself owns, all allocations that its ' .
        'inputs own, and allocations owned by software results this build uses.'
}
sub help_brief { return 'delays archiving of allocations used by these builds' };
sub help_synopsis { return help_brief() . "\n" };

sub execute {
    my $self = shift;
    my @allocations;
    for my $build ($self->builds) {
        push @allocations, $build->all_allocations;
    }

    Genome::Disk::Command::Allocation::DelayArchiving->execute(
        allocations => \@allocations,
        duration => $self->duration,
    );
    $self->status_message("Provided builds won't be archived for awhile!");
    return 1;
}

1;

