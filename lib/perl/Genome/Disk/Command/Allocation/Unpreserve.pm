package Genome::Disk::Command::Allocation::Unpreserve;

use strict;
use warnings;
use Genome;

class Genome::Disk::Command::Allocation::Unpreserve {
    is => 'Command::V2',
    has => [
        allocations => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'allocations to be unpreserved',
        },
    ],
    has_optional => [
        reason => {
            is => 'Text',
            doc => 'reason for wanting to set these allocations as unpreserved',
        },
    ],
    doc => 'unpreserves the given allocations',
};

sub help_detail {
    return 'unpreserves the given allocations, preventing them from being modified or deleted';
}
sub help_brief { return help_detail() }

sub execute {
    my $self = shift;
    $self->status_message("Starting de-preservation command!");

    for my $allocation ($self->allocations) {
        $self->debug_message("Unpreserving allocation " . $allocation->id);
        $allocation->preserved(0, $self->reason);
        unless (!$allocation->preserved) {
            Carp::confess "Could not unpreserve allocation " . $allocation->id;
        }
        $self->debug_message("Finished unpreserving allocation " . $allocation->id);
    }

    $self->status_message("Done unpreserving, exiting");
    return 1;
}

1;

