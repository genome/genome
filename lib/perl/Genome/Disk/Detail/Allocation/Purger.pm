package Genome::Disk::Detail::Allocation::Purger;

use strict;
use warnings;

use Genome;

use Carp qw(confess);
use File::Copy::Recursive qw(dirmove);

class Genome::Disk::Detail::Allocation::Purger {
    is => 'Genome::Disk::Detail::StrictObject',

    has => [
        allocation_id => {
            is => 'Text',
            len => 64,
        },
        reason => {
            is => 'Text',
        },
    ],

};

# Eventually SoftwareResults should have something like a status field,
# removing the need to set a test_name to invalidate a SoftwareResult
# backed by a purged allocation.
sub _update_owner_test_name {
    my $self       = shift;
    my $allocation = shift;
    my $event      = shift;

    my $owner = $allocation->owner;

    # Make sure the owner is a SoftwareResult before we set test_name.
    return unless defined $owner;
    return unless $owner->isa('Genome::SoftwareResult');
    return if defined $owner->test_name;

    my $timestamp = $event->created_at;
    my $reason    = $event->reason;

    my $test_name = sprintf "Allocation ID %s purged on %s, %s",
        $allocation->id, $timestamp,
        (defined $reason ? "with reason '$reason'" : 'no reason specified');

    $owner->test_name($test_name);

    return;
}

sub purge {
    my $self = shift;

    my $allocation_object = Genome::Disk::Allocation->get($self->allocation_id);
    die sprintf("No allocation found for id: %s",
        $self->allocation_id) unless $allocation_object;

    return 1 if $allocation_object->status eq 'purged';

    if ($allocation_object->is_archived) {
        return $self->_purge_archived($allocation_object);

    } else {
        return $self->_purge_unarchived($allocation_object);
    }
}


sub _purge_archived {
    my $self = shift;
    my $allocation_object = shift;

    unless ($ENV{UR_DBI_NO_COMMIT}) {
        my $archive_path = $allocation_object->archive_path;
        $self->status_message(q{Removing archive path '%s'}, $archive_path);
        unless (Genome::Sys->remove_directory_tree($archive_path)) {
            $self->error_message(
                'Error removing archive directory for allocation %s',
                $allocation_object->id
            );
            return;
        }
    }

    $self->_finalize_purge($allocation_object);

    return 1;
}


sub _purge_unarchived {
    my $self = shift;
    my $allocation_object = shift;

    unless ($ENV{UR_DBI_NO_COMMIT}) {
        my $absolute_path = $allocation_object->absolute_path;
        $self->status_message(q{Removing allocation path '%s'}, $absolute_path);
        unless (Genome::Sys->remove_directory_tree($absolute_path)) {
            $self->error_message(
                'Error removing directory for allocation %s',
                $allocation_object->id
            );
            return;
        }
    }

    $self->_finalize_purge($allocation_object);

    return 1;
}


sub _finalize_purge {
    my $self = shift;
    my $allocation_object = shift;

    my $event = Genome::Timeline::Event::Allocation->purged(
        $self->reason,
        $allocation_object,
    );

    $self->_update_owner_test_name($allocation_object, $event);

    $self->_update_allocation_status($allocation_object);

    return 1;
}


sub _update_allocation_status {
    my $self = shift;
    my $allocation_object = shift;

    $allocation_object->status('purged');
    $allocation_object->kilobytes_requested(0);
    $allocation_object->kilobytes_used(0);
}


1;
