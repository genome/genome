package Genome::Disk::Detail::Allocation::Archiver;

use strict;
use warnings;

use Genome;

use Carp qw(confess);
use Try::Tiny qw(try catch finally);

class Genome::Disk::Detail::Allocation::Archiver {
    is => 'Genome::Disk::Detail::StrictObject',

    has => [
        allocation_id => {
            is => 'Text',
            len => 64,
        },
    ],
};


sub archive {
    my $self = shift;
    my $id = $self->allocation_id;

    my ($allocation_object, $allocation_lock) =
        Genome::Disk::Allocation->get_with_lock($id);

    my $current_allocation_path = $allocation_object->absolute_path;
    my $archive_allocation_path = $allocation_object->archive_path();

    my $msg;
    my $creating_directory = 0;
    try {
        if ($allocation_object->is_archived) {
            confess sprintf("Allocation %s is already archived!",
                $allocation_object->id);
        }
        if (!$allocation_object->archivable) {
            confess sprintf("Allocation %s is not flagged as archivable!",
                $allocation_object->id);
        }
        unless (-e $current_allocation_path) {
            confess sprintf("Allocation path %s does not exist!",
                $current_allocation_path);
        }

        unless (Genome::Sys->recursively_validate_directory_for_read_write_access(
                $current_allocation_path)) {
            confess sprintf("Some files in allocation directory %s are not "
                . "readable or writable, this must be fixed before archiving "
                . "can continue", $current_allocation_path);
        }

        # Reallocate so we're reflecting the correct size at time of archive.
        $allocation_object->reallocate();

        unless ($allocation_object->volume->archive_is_mounted) {
            $self->fatal_message('volume for %s is not mounted', $allocation_object->archive_path);
        }

        unless ($allocation_object->volume->is_mounted) {
            $self->fatal_message('volume for %s is not mounted', $allocation_object->absolute_path);
        }

        $creating_directory = 1;
        Genome::Sys->create_directory($archive_allocation_path);
        Genome::Sys->rsync_directory(
            source_directory => $current_allocation_path,
            target_directory => $archive_allocation_path,
        );

        my $rv = $allocation_object->_commit_unless_testing;
        confess "Could not commit!" unless $rv;
    } catch {
        my $error = $_; # Record error so it can be investigated after unlocking

        my $cleanup_error;
        if ($creating_directory) {
            local $@;
            eval {
                Genome::Sys->remove_directory_tree($archive_allocation_path);
            };
            $cleanup_error = $@;
        }
        $msg = sprintf("Could not archive allocation %s, error:\n%s",
            $allocation_object->id, $error);
        if ($cleanup_error) {
            $msg .= sprintf("\n\nWhile cleaning up archive, "
                . "encoutered error:\n%s", $cleanup_error);
        }
    } finally {
        if ($allocation_lock) {
            $allocation_lock->unlock();
        }
    };
    confess $msg if $msg;

    Genome::Timeline::Event::Allocation->archived(
        'archived',
        $allocation_object,
    );
    $allocation_object->status('archived');

    # Never make filesystem changes if no commit is enabled
    unless ($ENV{UR_DBI_NO_COMMIT}) {
        my $rv = Genome::Sys->remove_directory_tree($current_allocation_path);
        unless ($rv) {
            confess sprintf("Could not remove unarchived allocation path %s. "
                . "Database changes have been committed, "
                . "so clean this up manually!", $current_allocation_path);
        }
    }

    return 1;
}


1;
