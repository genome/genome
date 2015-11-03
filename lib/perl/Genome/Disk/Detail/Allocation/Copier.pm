package Genome::Disk::Detail::Allocation::Copier;

use strict;
use warnings;

use Genome;

use Carp qw(confess);

class Genome::Disk::Detail::Allocation::Copier {
    is => 'Genome::Disk::Detail::StrictObject',

    has => [
        allocation_id => {
            is => 'Text',
            len => 64,
        },
        output_dir => {
            is => 'Text',
            len => 255,
        },
    ],
};

sub copy {
    my $self = shift;

    my $id = $self->allocation_id;

    my $output_dir = $self->output_dir;

    my ($allocation_object, $allocation_lock) =
        Genome::Disk::Allocation->get_with_lock($id);
    my $original_absolute_path = $allocation_object->absolute_path;
    
    # make shadow allocation
    my %creation_params = $self->_get_copy_shadow_params($allocation_object);

    # The shadow allocation is just a way of keeping track of our temporary
    # additional disk usage during the copy.
    my $shadow_allocation = Genome::Disk::Allocation->shadow_get_or_create(%creation_params);
    my $shadow_absolute_path = $shadow_allocation->absolute_path;
    my %rsync_params = (
        source_directory => $original_absolute_path,
        target_directory => $output_dir,
    );
    if ($allocation_object->is_archived) {
        Genome::Sys->create_directory($shadow_absolute_path);
        eval {
            my $tar_path = $allocation_object->tar_path;
            my $cmd = "tar -C $shadow_absolute_path -xf $tar_path";

            # It's very possible that if no commit is on, the volumes/allocations
            # being dealt with are test objects that don't exist out of this local
            # UR context, so bsubbing jobs would fail.
            if ($ENV{UR_DBI_NO_COMMIT}) {
                Genome::Sys->shellcmd(cmd => $cmd);
            }
            else {
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
                        die "Received signal, exiting.";
                    }
                }
                
                # Entire path must be wrapped in quotes because older allocation
                # IDs contain spaces If the command isn't wrapped in quotes, the
                # '&&' is misinterpreted by bash (rather than being "bsub '1 && 2'
                # it is looked at as 'bsub 1' && '2')
                ($job_id, $status) = Genome::Sys->bsub_and_wait(
                    queue => Genome::Config::get('archive_lsf_queue'),
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
            $rsync_params{source_directory} = $shadow_absolute_path;
        };
        my $unarchive_error_message = $@;
        if ($unarchive_error_message) {
            if ($shadow_absolute_path and -d $shadow_absolute_path and not $ENV{UR_DBI_NO_COMMIT}) {
                if (Genome::Sys->remove_directory_tree($shadow_absolute_path)) {
                    $shadow_allocation->delete;
                }
            }
            if ($allocation_lock) {
                $allocation_lock->unlock();
            }
            confess "Could not unarchive to shadow allocation, received error:\n$unarchive_error_message";
        }
    }

    # copy files to output_dir
    my $copy_rv = eval {
        Genome::Sys->rsync_directory(%rsync_params);
    };
    unless ($copy_rv) {
        my $copy_error_message = $@;
        if (-d $shadow_absolute_path) {
            if (Genome::Sys->remove_directory_tree($shadow_absolute_path)) {
                $shadow_allocation->delete;
            }
        }
        if ($allocation_lock) {
            $allocation_lock->unlock();
        }
        confess(sprintf(
            "Could not copy allocation %s from %s to %s: %s",
            $allocation_object->id, $original_absolute_path,
            $shadow_absolute_path, $copy_error_message));
    }

    Genome::Timeline::Event::Allocation->copied(
        sprintf("copied from %s to %s", $original_absolute_path,
                $output_dir),
        $allocation_object,
    );

    Genome::Disk::Allocation::_commit_unless_testing();

    $allocation_lock->unlock();

    $shadow_allocation->delete;

}

sub _get_copy_shadow_path {
    my $allocation_path = shift;
    return sprintf("%s-copy_allocation_destination", $allocation_path)
}

sub _get_copy_shadow_params {
    my ($self, $allocation) = @_;

    my %creation_parameters = (
        disk_group_name => $allocation->disk_group_name,
        kilobytes_requested => $allocation->kilobytes_requested,
        owner_class_name => "UR::Value",
        owner_id => "shadow_allocation",
        exclude_mount_path => $allocation->mount_path,
        allocation_path => _get_copy_shadow_path($allocation->allocation_path),
    );

    return %creation_parameters;
}


1;
