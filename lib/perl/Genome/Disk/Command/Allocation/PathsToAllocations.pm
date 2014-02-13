package Genome::Disk::Command::Allocation::PathsToAllocations;

use strict;
use warnings;

use Genome;

class Genome::Disk::Command::Allocation::PathsToAllocations {
    is => 'Command::V2',
    has => [
        absolute_paths => {
            is => 'DirectoryPath',
            is_many => 1,
            shell_args_position => 1,
            doc => 'Absolute path for which an allocation needs to be found',
        },
    ],
};

sub execute {
    my $self = shift;

    my @allocations;
    for my $path ($self->absolute_paths) {
        my $allocation = Genome::Disk::Allocation->get_allocation_for_path($path);
        if ($allocation) {
            push @allocations, $allocation;
        } else {
            $self->warning_message("Found no allocation for path $path");
        }
    }

    for my $allocation (@allocations) {
        print join("\t", $allocation->id, $allocation->absolute_path) . "\n";
    }
    return 1;
}

1;

