package Genome::Disk::Detail::Allocation::Archiver;

use strict;
use warnings;

use Genome;

use Carp qw(confess);

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
    my $archive_allocation_path = join('/',
        $allocation_object->volume->archive_mount_path,
        $allocation_object->group_subdirectory,
        $allocation_object->allocation_path);
    my $tar_path = join('/', $archive_allocation_path, 'archive.tar');

    # This gets set to true immediately before tarball creation is started.
    # This allows for conditional clean up of the archive directory in case of
    # failure, which is nice because that requires an LSF job be scheduled.
    my $tarball_created = 0;

    eval {
        if ($allocation_object->is_archived) {
            confess sprintf("Allocation %s is already archived!",
                $allocation_object->id);
        }
        if (!$allocation_object->archivable) {
            confess sprintf("Allocation %s is not falgged as archivable!",
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
        # Allocation size be greater than 1GB to be archived
        if (!$ENV{UR_DBI_NO_COMMIT}
                and $allocation_object->kilobytes_requested < 1048576) {
            confess(sprintf("Total size of files at path %s is only %s, "
                    . "which is not greater than 1GB!",
                    $current_allocation_path,
                    $allocation_object->kilobytes_requested));
        }

        my $mkdir_cmd = "mkdir -p $archive_allocation_path";
        my $cd_cmd = "cd $current_allocation_path";
        my $tar_cmd = "/bin/ls -A | tar --create --file $tar_path -T -";
        my $cmd = join(' && ', $mkdir_cmd, $cd_cmd, $tar_cmd);

        $tarball_created = 1;
        # It's very possible that if no commit is on, the volumes/allocations
        # being dealt with are test objects that don't exist out of this local
        # UR context, so bsubbing jobs would fail.
        if ($ENV{UR_DBI_NO_COMMIT}) {
            Genome::Sys->shellcmd(cmd => $cmd);
        }
        else {
            my ($job_id, $status);

            my @signals = qw/ INT TERM /;
            for my $signal (@signals) {
                $SIG{$signal} = sub {
                    print STDERR "Cleanup activated within allocation, cleaning up LSF jobs\n";
                    eval { Genome::Sys->kill_lsf_job($job_id) } if $job_id;
                    $allocation_object->_cleanup_archive_directory($archive_allocation_path);
                    die "Received signal, exiting.";
                }
            }
            # Entire path must be wrapped in quotes because older allocation
            # IDs contain spaces If the command isn't wrapped in quotes, the
            # '&&' is misinterpreted by bash (rather than being "bsub '1 && 2'
            # it is looked at as 'bsub 1' && '2')
            ($job_id, $status) = Genome::Sys->bsub_and_wait(
                queue => $ENV{GENOME_ARCHIVE_LSF_QUEUE},
                job_group => '/archive',
                log_file => "/tmp/$id",
                cmd => "$cmd",
            );

            for my $signal (@signals) {
                delete $SIG{$signal};
            }

            unless ($status eq 'DONE') {
                confess sprintf("LSF job %s failed to execute %s, "
                    . "exited with status %s", $job_id, $cmd, $status);
            }
        }

        $allocation_object->mount_path(
            $allocation_object->volume->archive_mount_path);
        $allocation_object->_update_owner_for_move;

        my $rv = $allocation_object->_commit_unless_testing;
        confess "Could not commit!" unless $rv;
    };
    my $error = $@; # Record error so it can be investigated after unlocking

    # If only there were finally blocks...
    if ($allocation_lock) {
        Genome::Sys->unlock_resource(resource_lock => $allocation_lock);
    }

    if ($error) {
        eval {
            $allocation_object->_cleanup_archive_directory(
                $archive_allocation_path) if $tarball_created;
        };
        my $cleanup_error = $@;
        my $msg = sprintf("Could not archive allocation %s, error:\n%s",
            $allocation_object->id, $error);
        if ($cleanup_error) {
            $msg .= sprintf("\n\nWhile cleaning up archive tarball, "
                . "encoutered error:\n%s", $cleanup_error);
        }
        confess $msg;
    } else {
        Genome::Timeline::Event::Allocation->archived(
            'archived',
            $allocation_object,
        );
        $allocation_object->status('archived');
    }

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
