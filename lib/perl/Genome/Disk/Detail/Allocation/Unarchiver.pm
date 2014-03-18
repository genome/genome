package Genome::Disk::Detail::Allocation::Unarchiver;

use strict;
use warnings;

use Genome;

use Carp qw(confess);

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
        Genome::Sys->unlock_resource(resource_lock => $allocation_lock);
        $allocation_object->status_message(
            "Allocation is not archived, cannot unarchive. Exiting.");
        return 1;
    }

    my %creation_params = $allocation_object->unarchive_shadow_params;
    # shadow_allocation ensures that we wont over allocate our destination volume
    my $shadow_allocation = Genome::Disk::Allocation->shadow_get_or_create(
        %creation_params);

    my $archive_path = $allocation_object->absolute_path;
    my $target_path = $shadow_allocation->absolute_path;

    # Inferred path prior to archiving so we can symlink the new allocation,
    # to this old location.
    (my $old_absolute_path = $archive_path) =~ s/^\/gscarchive/\/gscmnt/;

    my $tar_path = $allocation_object->tar_path;
    my $cmd = "tar -C $target_path -xf $tar_path";

    eval {
        # It's very possible that if no commit is on, the volumes/allocations
        # being dealt with are test objects that don't exist out of this local
        # UR context, so bsubbing jobs would fail.
        if ($ENV{UR_DBI_NO_COMMIT}) {
            Genome::Sys->shellcmd(cmd => $cmd);
        } else {
            # If this process should be killed, the LSF job needs to be cleaned up
            my ($job_id, $status);

            # Signal handlers are added like this so an anonymous sub can be
            # used, which handles variables defined in an outer scope
            # differently than named subs (in this case, $job_id,
            # $allocation_object, and $archive_path).
            my @signals = qw/ INT TERM /;
            for my $signal (@signals) {
                $SIG{$signal} = sub {
                    print STDERR "Cleanup activated within allocation, ",
                        "cleaning up LSF jobs\n";
                    eval { Genome::Sys->kill_lsf_job($job_id) } if $job_id;
                    $allocation_object->_cleanup_archive_directory(
                        $archive_path);
                    die "Received signal, exiting.";
                }
            }

            # Entire path must be wrapped in quotes because older allocation
            # IDs contain spaces If the command isn't wrapped in quotes, the
            # '&&' is misinterpreted by bash (rather than being "bsub '1 && 2'
            # it is looked at as 'bsub 1' && '2')
            ($job_id, $status) = Genome::Sys->bsub_and_wait(
                queue => $ENV{GENOME_ARCHIVE_LSF_QUEUE},
                job_group => '/unarchive',
                log_file => "/tmp/$id",
                cmd => $cmd,
            );

            for my $signal (@signals) {
                delete $SIG{$signal};
            }

            unless ($status eq 'DONE') {
                confess "Could not execute command $cmd via LSF job $job_id, "
                    . "received status $status";
            }
        }

        # Make updates to the allocation
        $allocation_object->mount_path($shadow_allocation->volume->mount_path);
        Genome::Sys->create_directory($allocation_object->absolute_path);
        unless (rename $shadow_allocation->absolute_path,
                $allocation_object->absolute_path) {
            confess($allocation_object->error_message(sprintf(
                    "Could not move shadow allocation path (%s) to "
                    . "final path (%s).  This should never happen, "
                    . "even when 100%% full.",
                    $shadow_allocation->absolute_path,
                    $allocation_object->absolute_path)));
        }
        $allocation_object->_update_owner_for_move;
        $allocation_object->archive_after_time(
            Genome::Disk::Command::Allocation::DelayArchiving->_resolve_date_from_months(3));

        if ($old_absolute_path ne $allocation_object->absolute_path) {
            Genome::Disk::Allocation::_symlink_new_path_from_old(
                $old_absolute_path, $allocation_object->absolute_path);
        }

        unless ($allocation_object->_commit_unless_testing) {
            confess "Could not commit!";
        }
    };
    my $error = $@;

    # finally blocks would be really sweet. Alas...
    Genome::Sys->unlock_resource(resource_lock => $allocation_lock) if $allocation_lock;
    $shadow_allocation->delete();

    if ($error) {
        if ($target_path and -d $target_path and not $ENV{UR_DBI_NO_COMMIT}) {
            Genome::Sys->remove_directory_tree($target_path);
        }
        confess "Could not unarchive, received error:\n$error";
    } else {
        $allocation_object->add_note(
            header_text => 'unarchived',
            body_text => $self->reason,
        );
        Genome::Timeline::Event::Allocation->unarchived($self->reason,
            $allocation_object);
        $allocation_object->status('active');
    }

    $allocation_object->_cleanup_archive_directory($archive_path);
    return 1;
}

1;
