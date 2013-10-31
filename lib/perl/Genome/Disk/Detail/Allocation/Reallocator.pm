package Genome::Disk::Detail::Allocation::Reallocator;

use strict;
use warnings;

use Genome;

use Carp qw(confess);

class Genome::Disk::Detail::Allocation::Reallocator {
    is => 'Genome::Disk::Detail::StrictObject',

    has => [
        allocation_id => {
            is => 'Text',
            len => 64,
        },
    ],

    has_optional => [
        kilobytes_requested => {
            is => 'Number',
            len => 20,
        },

        allow_reallocate_with_move => {
            is => 'Boolean',
            default => 0,
        },

        grow_only => {
            is => 'Boolean',
            default => 0,
        },
    ],
};


sub reallocate {
    my $self = shift;

    my $mode = Genome::Disk::Allocation->_retrieve_mode();
    my $allocation_object = Genome::Disk::Allocation->$mode(
        $self->allocation_id);
    my $old_kb_requested = $allocation_object->kilobytes_requested;

    my $kb_used = $allocation_object->du();

    # Cache kilobytes used
    $allocation_object->kilobytes_used($kb_used);

    my $actual_kb_requested = List::Util::max($kb_used,
        $self->kilobytes_requested);
    if ($self->grow_only && ($actual_kb_requested <= $old_kb_requested)) {
        $allocation_object->status_message(
            "Not changing kilobytes_requested, because grow_only = 1 & actual usage < original_kilobytes_requested");
        return 1;
    }
    if ($actual_kb_requested > $self->kilobytes_requested) {
        $allocation_object->status_message(sprintf(
                "Setting kilobytes_requested to %s based on `du` for allocation %s",
                $actual_kb_requested, $allocation_object->id));
    }

    $allocation_object->reallocation_time(UR::Context->current->now);
    $allocation_object->kilobytes_requested($actual_kb_requested);
    Genome::Disk::Allocation::_commit_unless_testing();

    my $volume = $allocation_object->volume;
    my $succeeded;
    if ($volume->allocated_kb < $volume->hard_limit_kb) {
        $succeeded = 1;
    } else {
        if ($self->allow_reallocate_with_move) {
            my $old_mount_path = $allocation_object->mount_path;
            $allocation_object->move();
            unless ($allocation_object->mount_path eq $old_mount_path) {
                $succeeded = 1;
            } else {
                $allocation_object->error_message(sprintf(
                        "Failed to reallocate and move allocation %s from volume %s.",
                        $allocation_object->id, $volume->mount_path));
            }
        } else {
            $allocation_object->error_message(sprintf(
                    "Failed to reallocate allocation %s on volume %s because the volume is beyond its quota.",
                    $allocation_object->id, $volume->mount_path));
        }
    }

    if ($succeeded) {
        Genome::Timeline::Event::Allocation->reallocated(
            "actual kb used: $kb_used",
            $allocation_object,
        );
        Genome::Disk::Allocation::_commit_unless_testing();
    } else {
        # Rollback kilobytes_requested
        my $max_kilobytes_requested = List::Util::max($kb_used, $old_kb_requested);
        my $msg = $old_kb_requested == $max_kilobytes_requested ? 'Rolling back' : 'Setting';

        $allocation_object->status_message(sprintf(
                "%s kilobytes_requested to %d for allocation %s.",
                $msg, $max_kilobytes_requested, $allocation_object->id));

        $allocation_object->kilobytes_requested($max_kilobytes_requested);
        Genome::Disk::Allocation::_commit_unless_testing();
    }

    return $succeeded;

}


1;
