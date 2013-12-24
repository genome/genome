package Genome::Disk::Detail::Allocation::CreationParameters;

use strict;
use warnings;

use Genome;
use UR;

use Carp qw(confess);

class Genome::Disk::Detail::Allocation::CreationParameters {
    is => 'Genome::Disk::Detail::StrictObject',
    has => [
        kilobytes_requested => {
            is => 'Number',
            len => 20,
        },

        owner_class_name => {
            is => 'Text',
            len => 255,
        },
        owner_id => {
            is => 'Text',
            len => 255,
        },

        allocation_path => {
            is => 'Text',
            len => 4000,
        },
        disk_group_name => {
            is => 'Text',
            len => 40,
        },
    ],

    has_optional => [
        allocation_id => {
            is => 'Text',
            len => 64,
        },

        group_subdirectory => {
            is => 'Text',
            len => 255,
        },

        mount_path => {
            is => 'Text',
            len => 255,
        },
        exclude_mount_path => {
            is => 'Text',
            len => 255,
        },

        archive_after_time => {
            is => 'DateTime',
            len => 11,
        },
        kilobytes_used => {
            is => 'Number',
            len => 20,
        },
    ],
};


sub get_id {
    my $self = shift;
    if ($self->allocation_id) {
        return $self->allocation_id;
    } else {
        # TODO autogenerate_new_object_id should technically receive a BoolExpr
        return Genome::Disk::Allocation->__meta__->autogenerate_new_object_id;
    }
}

sub sanitize {
    my $self = shift;

    $self->group_subdirectory($self->disk_group->subdirectory);
    $self->sanitize_mount_path;
    $self->sanitize_exclude_mount_path;
}

sub sanitize_mount_path {
    my $self = shift;

    $self->mount_path(sanitize_directory_path($self->mount_path));
}

sub sanitize_exclude_mount_path {
    my $self = shift;

    $self->exclude_mount_path(sanitize_directory_path(
            $self->exclude_mount_path));
}

sub sanitize_directory_path {
    my $path = shift;
    if ($path) {
        $path =~ s/\/$//; # mount paths in database don't have trailing /
    }
    return $path;
}

sub validate {
    my $self = shift;

    $self->validate_owner_class_name;
    $self->validate_kilobytes_requested;
    $self->validate_group_subdirectory;
    $self->validate_mount_path;
}

sub validate_owner_class_name {
    my $self = shift;

    unless ($self->owner_class_name->__meta__) {
        confess sprintf("Could not find meta information for owner class %s, "
            . "make sure this class exists!", $self->owner_class_name);
    }
}

sub validate_kilobytes_requested {
    my $self = shift;

    unless ($self->kilobytes_requested >= 0) {
        confess sprintf('Kilobytes requested is negative (%s)!',
            $self->kilobytes_requested);
    }
}

sub validate_group_subdirectory {
    my $self = shift;

    my $group = $self->disk_group;

    if (defined $self->group_subdirectory
            and $self->group_subdirectory ne $group->subdirectory) {
        $self->warning_message(sprintf(
            "Given group subdirectory %s does not match retrieved "
            . "group's subdirectory, ignoring provided value\n",
            $self->group_subdirectory));
    }
}

sub validate_mount_path {
    my $self = shift;
    if (defined $self->mount_path && defined $self->exclude_mount_path) {
        if ($self->mount_path eq $self->exclude_mount_path) {
            confess sprintf("mount_path (%s) equal to exclude_mount_path (%s)",
                $self->mount_path, $self->exclude_mount_path);
        }
    }
}

sub disk_group {
    my $self = shift;

    my $group = Genome::Disk::Group->get(
        disk_group_name => $self->disk_group_name);
    confess sprintf("Could not find a group with name %s",
        $self->disk_group_name) unless $group;
    return $group;
}

sub as_hash {
    my $self = shift;

    my %parameters = (
        disk_group_name              => $self->disk_group_name,
        kilobytes_requested          => $self->kilobytes_requested,
        original_kilobytes_requested => $self->kilobytes_requested,
        allocation_path              => $self->allocation_path,
        owner_class_name             => $self->owner_class_name,
        owner_id                     => $self->owner_id,
        group_subdirectory           => $self->group_subdirectory,
        creation_time                => UR::Context->current->now,
    );
    if ($self->allocation_id) {
        $parameters{id} = $self->allocation_id;
    }
    if ($self->archive_after_time) {
        $parameters{archive_after_time} = $self->archive_after_time;
    }
    if ($self->kilobytes_used) {
        $parameters{kilobytes_used} = $self->kilobytes_used;
    }
    return %parameters;
}

1;
