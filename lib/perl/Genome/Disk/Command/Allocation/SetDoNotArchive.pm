package Genome::Disk::Command::Allocation::SetDoNotArchive;

use strict;
use warnings;
use Genome;

class Genome::Disk::Command::Allocation::SetDoNotArchive {
    is => 'Command::V2',
    has_optional => [
        allocations => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            doc => 'allocations that are not to be archived, resolved via Command::V2',
        },
        paths => {
            is => 'Text',
            doc => 'comma delimited list of paths, an attempt will be made to resolve them to allocations',
        },
        reason => {
            is => 'Text',
            doc => 'reason for wanting these allocations/paths to not be archived',
        },
    ],
};

sub help_detail { 
    return 'Given allocations and paths are set so they cannot be archived. Any paths that are given are resolved ' .
        'to allocations if possible, otherwise a warning is emitted.';
}
sub help_brief { return 'given allocations and paths are marked so they cannot be archived' };
sub help_synopsis { return help_brief() . "\n" };

sub execute {
    my $self = shift;
    for my $allocation ($self->_resolve_allocations_from_paths, $self->allocations) {
        next if !$allocation->archivable;
        $allocation->archivable(0, $self->reason);
    }
    $self->status_message("Allocations can no longer be archived");
    return 1;
}

# TODO Move this logic into a special get method on Genome::Disk::Allocation?
sub _resolve_allocations_from_paths {
    my $self = shift;
    return unless $self->paths;

    my @allocations;
    for my $path (split(/,/, $self->paths)) {
        my $allocation_path = Genome::Disk::Allocation->_allocation_path_from_full_path($path);
        unless ($allocation_path) {
            $self->warning_message("No allocation found for path $path");
            next;
        }

        my $allocation = Genome::Disk::Allocation->get_parent_allocation($allocation_path);
        if ($allocation) {
            push @allocations, $allocation;
            next;
        }

        my @allocations_for_path = Genome::Disk::Allocation->get_child_allocations($allocation_path);
        if (@allocations_for_path) {
            push @allocations, @allocations_for_path;
            next;
        }

        $self->warning_message("No allocation found for path $path");
    }
    return @allocations;
}

1;

