package Genome::Disk::Detail::Allocation::Creator;

use strict;
use warnings;

use Genome;

use Carp qw(croak confess);


class Genome::Disk::Detail::Allocation::Creator {
    has => {
        parameters => {
            is => 'Genome::Disk::Detail::Allocation::CreationParameters',
        },
    },
};


sub create_allocation {
    my $self = shift;

    my $class = 'Genome::Disk::Allocation';

    # Make sure there aren't any extra params
    my $id = $self->parameters->get_id;


    my $kilobytes_requested = $self->parameters->kilobytes_requested;
    my $owner_class_name = $self->parameters->owner_class_name;
    my $owner_id = $self->parameters->owner_id;
    my $allocation_path = $self->parameters->allocation_path;
    my $disk_group_name = $self->parameters->disk_group_name;
    my $mount_path = $self->parameters->mount_path;
    my $exclude_mount_path = $self->parameters->exclude_mount_path;
    my $group_subdirectory = $self->parameters->group_subdirectory;

    if (my $parent_alloc = $class->get_parent_allocation($allocation_path)) {
        confess sprintf("Parent allocation (%s) found for %s", $parent_alloc->allocation_path, $allocation_path);
    }
    unless ($class->_verify_no_child_allocations($allocation_path)) {
        confess "Child allocation found for $allocation_path!";
    }

    if ($ENV{GENOME_DB_PAUSE} and -e $ENV{GENOME_DB_PAUSE}) {
        print "Database updating has been paused; not going to attempt to allocate disk until the pause is released. Please stand by...\n";

        while (1) {
            sleep 30;
            last unless -e $ENV{GENOME_DB_PAUSE};
        }

        print "Database updating has been resumed, continuing allocation!\n";
    }

    my $group = $self->parameters->disk_group;
    if (defined $group_subdirectory and $group_subdirectory ne $group->subdirectory) {
        print STDERR "Given group subdirectory $group_subdirectory does not match retrieved group's subdirectory, ignoring provided value\n";
    }
    $group_subdirectory = $group->subdirectory;

    # If given a mount path, need to ensure it's valid by trying to get a disk volume with it. Also need to make
    # sure that the retrieved volume actually belongs to the supplied disk group and that it can be allocated to
    my @candidate_volumes;
    Genome::Utility::Instrumentation::timer('disk.allocation.create.candidate_volumes.selection', sub {
        if (defined $mount_path) {
            $mount_path =~ s/\/$//; # mount paths in database don't have trailing /
            my $volume = Genome::Disk::Volume->get(mount_path => $mount_path, disk_status => 'active', can_allocate => 1);
            confess "Could not get volume with mount path $mount_path" unless $volume;

            unless (grep { $_ eq $disk_group_name } $volume->disk_group_names) {
                confess "Volume with mount path $mount_path is not in supplied group $disk_group_name!";
            }

            my @reasons;
            push @reasons, 'disk is not active' if $volume->disk_status ne 'active';
            push @reasons, 'allocation turned off for this disk' if $volume->can_allocate != 1;

            if ($exclude_mount_path && $volume->mount_path eq $exclude_mount_path) {
                push @reasons, 'Specified mount path matched the excluded mount path.';
            }

            if (@reasons) {
                confess "Requested volume with mount path $mount_path cannot be allocated to:\n" . join("\n", @reasons);
            }

            push @candidate_volumes, $volume;
        } else {
            my %candidate_volume_params = (
                disk_group_name => $disk_group_name,
            );
            if (defined $exclude_mount_path) {
                $candidate_volume_params{'exclude'} = $exclude_mount_path;
            }
            push @candidate_volumes, $class->_get_candidate_volumes(
                %candidate_volume_params);
        }
    });

    my %parameters = (
        disk_group_name              => $disk_group_name,
        kilobytes_requested          => $kilobytes_requested,
        original_kilobytes_requested => $kilobytes_requested,
        allocation_path              => $allocation_path,
        owner_class_name             => $owner_class_name,
        owner_id                     => $owner_id,
        group_subdirectory           => $group_subdirectory,
        id                           => $id,
        creation_time                => UR::Context->current->now,
    );
    if ($self->parameters->archive_after_time) {
        $parameters{archive_after_time} = $self->parameters->archive_after_time;
    }
    if ($self->parameters->kilobytes_used) {
        $parameters{kilobytes_used} = $self->parameters->kilobytes_used;
    }

    my $allocation_object;
    Genome::Utility::Instrumentation::timer('disk.allocation.create.get_allocation_without_lock', sub {
        $allocation_object = $class->_get_allocation_without_lock(\@candidate_volumes, \%parameters);
    });

    $allocation_object->debug_message(sprintf("Allocation (%s) created at %s",
        $id, $allocation_object->absolute_path));

    # a restrictive umask can break builds for other users; force it to be friendly
    umask(0002);
    # If we cannot create the directory delete the new allocation
    my $dir;
    eval {
        Genome::Utility::Instrumentation::timer('disk.allocation.create.create_directory', sub {
            $dir = Genome::Sys->create_directory($allocation_object->absolute_path);
        });
    };
    my $error = $@;
    unless (defined($dir) and ( -d $dir ) and not $error) {
        $class->error_message(sprintf(
                "Failed to create directory (%s) with return value = '%s', and error:\n%s",
                $allocation_object->absolute_path, $dir || '', $error));
        $allocation_object->delete;
        confess $error;
    }

    Genome::Timeline::Event::Allocation->created('initial creation', $allocation_object);

    return $allocation_object;
}

1;
