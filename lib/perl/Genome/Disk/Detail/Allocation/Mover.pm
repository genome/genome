package Genome::Disk::Detail::Allocation::Mover;

use strict;
use warnings;

use Genome;

use Carp qw(confess);

class Genome::Disk::Detail::Allocation::Mover {
    is => 'Genome::Disk::Detail::StrictObject',

    has => [
        allocation_id => {
            is => 'Text',
            len => 64,
        },
    ],

    has_optional => [
        disk_group_name => {
            is => 'Text',
            len => 40,
        },

        target_mount_path => {
            is => 'Text',
            len => 255,
        },
    ],
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    if ( $self->disk_group_name and $self->target_mount_path ) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => [qw/ group_name target_mount_path /],
            desc => 'Can not specify both target group and volume!', 
        );
    }

    return @errors;
}

sub move {
    my $self = shift;

    my $id = $self->allocation_id;

    my $group_name = $self->disk_group_name;
    my $new_mount_path = $self->target_mount_path;

    my ($allocation_object, $allocation_lock) =
        Genome::Disk::Allocation->get_with_lock($id);

    my $original_absolute_path = $allocation_object->absolute_path;

    # make shadow allocation
    my %creation_params = $self->_get_move_shadow_params($allocation_object);

    # The shadow allocation is just a way of keeping track of our temporary
    # additional disk usage during the move.
    my $shadow_allocation = Genome::Disk::Allocation->shadow_get_or_create(
        %creation_params);

    my $shadow_absolute_path = $shadow_allocation->absolute_path;

    # copy files to shadow allocation
    my $copy_rv = eval {
        Genome::Sys->rsync_directory(
            source_directory => $original_absolute_path,
            target_directory => $shadow_absolute_path,
        );
    };
    unless ($copy_rv) {
        my $copy_error_message = $@;
        if (-d $shadow_absolute_path) {
            if (Genome::Sys->remove_directory_tree($shadow_absolute_path)) {
                $shadow_allocation->delete;
            }
        }
        $allocation_lock->unlock();
        confess(sprintf(
                "Could not copy allocation %s from %s to %s: %s",
                $allocation_object->id, $original_absolute_path,
                $shadow_absolute_path, $copy_error_message));
    }

    my $new_volume_final_path = Genome::Disk::Allocation->_absolute_path(
        $shadow_allocation->mount_path,
        $shadow_allocation->group_subdirectory,
        $allocation_object->allocation_path
    );

    if (-e $new_volume_final_path) {
        if (-d $new_volume_final_path) {
            # Recursively compare directory paths, only collecting different file names
            my $diff = `diff -rq $original_absolute_path $new_volume_final_path`;
            if ($diff) {
                $allocation_lock->unlock();
                $shadow_allocation->delete;
                confess(
                    $allocation_object->error_message(
                        sprintf(
                            "Could not move original allocation path (%s) to final path (%s)."
                                . "  Differences found between directories: %s",
                            $shadow_absolute_path,
                            $new_volume_final_path,
                            $diff,
                        )
                    )
                );
            }
            $self->warning_message(
                "removing existing path (%s), confirmed redundant with original path (%s)",
                $new_volume_final_path,
                $original_absolute_path
            );
            # No differences detected, remove new path to allow rename of shadow allocation
            unless (Genome::Sys->remove_directory_tree($new_volume_final_path)) {
                $allocation_lock->unlock();
                $shadow_allocation->delete;
                confess(
                    $allocation_object->error_message(
                        sprintf(
                            "Could not remove existing, redundant final path (%s).",
                            $new_volume_final_path
                        )
                    )
                );
            }
        } else {
            # Something exists, but not a directory
            $allocation_lock->unlock();
            $shadow_allocation->delete;
            confess(
                $allocation_object->error_message(
                    sprintf(
                        "Unexpected non-directory path (%s) already exists!",
                        $new_volume_final_path
                    )
                )
            );
        }
    }

    Genome::Sys->create_directory($new_volume_final_path);
    unless (Genome::Sys->rename($shadow_allocation->absolute_path, $new_volume_final_path)) {
        $allocation_lock->unlock();
        $shadow_allocation->delete;
        confess(
            $allocation_object->error_message(
                sprintf(
                    "Could not move shadow allocation path (%s) to final path (%s)."
                        . "  This should never happen, even when 100%% full.",
                    $shadow_absolute_path,
                    $new_volume_final_path
                )
            )
        );
    }

    # Change the shadow allocation to reserve some disk on the old volume until
    # those files are deleted.
    my $old_mount_path = $allocation_object->mount_path;
    $allocation_object->mount_path($shadow_allocation->mount_path);
    $allocation_object->group_subdirectory($shadow_allocation->group_subdirectory);
    $allocation_object->disk_group_name($shadow_allocation->disk_group_name);
    $shadow_allocation->mount_path($old_mount_path);

    # This is here because certain objects (Build & SoftwareResult) don't
    # calculate their data_directories from their disk_allocations.
    $allocation_object->_update_owner_for_move;

    Genome::Timeline::Event::Allocation->moved(
        sprintf("moved from %s to %s", $original_absolute_path,
            $allocation_object->absolute_path),
        $allocation_object,
    );

    Genome::Disk::Allocation::_commit_unless_testing();

    $allocation_lock->unlock();

    $shadow_allocation->delete;

    Genome::Disk::Allocation->_create_observer(
        Genome::Disk::Allocation->_mark_for_deletion_closure(
            $original_absolute_path),
        Genome::Disk::Allocation->_remove_directory_closure(
            $original_absolute_path),
        sub { Genome::Disk::Allocation::_symlink_new_path_from_old(
                $allocation_object->absolute_path, $original_absolute_path) },
    );
}

sub _get_move_shadow_path {
    my $allocation_path = shift;
    return sprintf("%s-move_allocation_destination", $allocation_path)
}

sub _get_move_shadow_params {
    my ($self, $allocation) = @_;

    my %creation_parameters = (
        disk_group_name => $allocation->disk_group_name,
        kilobytes_requested => $allocation->kilobytes_requested,
        owner_class_name => "UR::Value",
        owner_id => "shadow_allocation",
        exclude_mount_path => $allocation->mount_path,
        allocation_path => _get_move_shadow_path($allocation->allocation_path),
    );

    # I think that it's dangerous to specify the new mount path, but this
    # feature existed, so nnutter and I kept it during this refactor.
    if ($self->target_mount_path) {
        $creation_parameters{'mount_path'} = $self->target_mount_path;
    }

    if ($self->disk_group_name) {
        $creation_parameters{disk_group_name} = $self->disk_group_name;
    }


    return %creation_parameters;
}


1;
