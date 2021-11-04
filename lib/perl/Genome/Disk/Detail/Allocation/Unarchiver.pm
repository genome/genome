package Genome::Disk::Detail::Allocation::Unarchiver;

use strict;
use warnings;

use Genome;

use Carp qw(confess);
use Try::Tiny;

class Genome::Disk::Detail::Allocation::Unarchiver {
    is => 'Genome::Disk::Detail::StrictObject',

    has => [
        allocation_id => {
            is => 'Text',
            len => 64,
        },
    ],

    has_optional => [
        reason => {
            is => 'Text',
            default_value => 'no reason given',
        },
    ],
};


sub unarchive {
    my $self = shift;

    my $id = $self->allocation_id;

    my ($allocation_object, $allocation_lock) =
        Genome::Disk::Allocation->get_with_lock($id);

    unless ($allocation_object->is_archived) {
        $allocation_lock->unlock();
        $allocation_object->status_message(
            "Allocation is not archived, cannot unarchive. Exiting.");
        return 1;
    }

    if ($allocation_object->mount_path =~ m!^/gscarchive/!) {
        $allocation_lock->unlock();
        $allocation_object->error_message('Cannot unarchive data from old archive system.');
        return;
    }

    my $tx = UR::Context::Transaction->begin();
    try {
        $self->_do_unarchive_cmd($allocation_object);

        $allocation_object->archive_after_time(
            Genome::Disk::Command::Allocation::DelayArchiving->_resolve_date_from_months(3));

        $allocation_object->add_note(
            header_text => 'unarchived',
            body_text => $self->reason,
        );
        Genome::Timeline::Event::Allocation->unarchived($self->reason, $allocation_object);
        $allocation_object->status('active');

        if ($tx->commit()) {
            undef $tx;
            unless ($allocation_object->_commit_unless_testing) {
                die 'failed to commit';
            }
        }
        else {
            die 'failed to commit transaction';
        }
    }
    catch {
        my $error = $_;
        if ($tx) {
            $tx->rollback();
        }
        my $target_path = $allocation_object->absolute_path;
        if ($target_path and -d $target_path and not $ENV{UR_DBI_NO_COMMIT}) {
            Genome::Sys->remove_directory_tree($target_path);
        }
        $self->fatal_message('Could not unarchive, received error: %s', $error);
    }
    finally {
        $allocation_lock->unlock() if $allocation_lock;
    };

    unless ($ENV{UR_DBI_NO_COMMIT}) {
        my $rv = Genome::Sys->remove_directory_tree($allocation_object->archive_path);
        unless ($rv) {
            $self->fatal_message('Could not remove archive for allocation %s; manual cleanup may be required.', $allocation_object->id);
        }
    }

    return 1;
}

sub _do_unarchive_cmd {
    my $class = shift;
    my $archived_allocation = shift;

    my $active_allocation_path = $archived_allocation->absolute_path;
    my $archive_allocation_path = $archived_allocation->archive_path();

    unless ($archived_allocation->volume->archive_is_mounted) {
        $class->fatal_message('volume for %s is not mounted', $archived_allocation->archive_path);
    }

    unless ($archived_allocation->volume->is_mounted) {
        $class->fatal_message('volume for %s is not mounted', $archived_allocation->absolute_path);
    }

    Genome::Sys->create_directory($active_allocation_path);
    Genome::Sys->rsync_directory(
        source_directory => $archive_allocation_path,
        target_directory => $active_allocation_path,
    );

    return 1;
}

1;
