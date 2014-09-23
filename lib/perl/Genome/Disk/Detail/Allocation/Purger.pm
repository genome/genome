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
        $allocation_object->_cleanup_archive_directory(
            $allocation_object->absolute_path);
    }

    $self->_finalize_purge($allocation_object);

    return 1;
}


sub _purge_unarchived {
    my $self = shift;
    my $allocation_object = shift;

    $allocation_object->_create_file_summaries();

    unless ($ENV{UR_DBI_NO_COMMIT}) {
        my $destination_directory = $self->_destination_dir($allocation_object);
        $self->status_message('Moving allocation path \''.
            $allocation_object->absolute_path .'\' to temporary path \''.
            $destination_directory .'\'');

        unless (dirmove($allocation_object->absolute_path,
                    $destination_directory)) {
            $self->error_message('Failed to move allocation path \''.
                $allocation_object->absolute_path .'\' to destination path \''.
                $destination_directory .'\': '. $!);
            return;
        };
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

    $allocation_object->status('purged');
    $allocation_object->kilobytes_requested(0);
    $allocation_object->kilobytes_used(0);

    return 1;
}


sub _destination_dir {
    my $self = shift;
    my $allocation_object = shift;

    return Genome::Sys->create_directory(File::Spec->join(
            $allocation_object->_get_trash_folder(),
            $allocation_object->id));
}


1;
