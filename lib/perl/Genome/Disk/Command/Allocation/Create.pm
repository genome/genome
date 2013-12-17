package Genome::Disk::Command::Allocation::Create;

use strict;
use warnings;
use Genome;

class Genome::Disk::Command::Allocation::Create {
    is => 'Command::V2',
    has => [
        disk_group_name => {
            is => 'Text',
            doc => 'Name of the disk group in which the allocation should be made.',
            valid_values => [$ENV{GENOME_DISK_GROUP_DEV}, $ENV{GENOME_DISK_GROUP_REFERENCES}, $ENV{GENOME_DISK_GROUP_ALIGNMENTS}, $ENV{GENOME_DISK_GROUP_MODELS}, $ENV{GENOME_DISK_GROUP_RESEARCH}],
        },
        allocation_path => {
            is => 'Text',
            doc => 'Subdirectory of the disk volume for which space is being allocated.',
        },
        kilobytes_requested => {
            is => 'Number',
            doc => 'Number of kilobytes to be reserved.',
        },
        owner_class_name => {
            is => 'Text',
            default_value => 'Genome::Sys::User',
            doc => 'Class name of entity that owns the allocations (eg, Genome::Sys::User for users, or Genome::Model::Build::* for builds.',
        },
        owner_id => {
            is => 'Text',
            default_value => Genome::Sys::User->owner_id,
            doc => 'The ID used to retrieve the owner (in conjunction with owner_class_name), e.g. username@example.com for Genome::Sys::User or build_id for Genome::Model::Build.',
        },
    ],
    has_optional => [
        mount_path => {
            is => 'Text',
            doc => 'The mount path for the disk volume on which the allocation should be created, if none is given then one is found and assigned automatically.',
        },
        kilobytes_used => {
            is => 'Number',
            doc => 'The amount of space currently consumed by data within the allocated directory, updated by reallocation.',
            default => 0,
        },
    ],
    doc => 'Creates an allocation for data somewhere on the filesystem.'
};

sub help_brief {
    return 'creates an allocation for data somewhere on the filesystem';
}

sub help_synopsis {
    return 'Creates an allocation for data somewhere on the filesystem, which reserves the space for the owner';
}

sub help_detail {
    return <<EOS
This tool creates an allocation for space somewhere on the filesystem. The user can specify
exactly where by providing the mount path and allocation path, or a mount path can be left
out and picked automatically.
EOS
}

sub execute {
    my $self = shift;
    my %params = (
        disk_group_name => $self->disk_group_name,
        allocation_path => $self->allocation_path,
        kilobytes_requested => $self->kilobytes_requested,
        owner_class_name => $self->owner_class_name,
        owner_id => $self->owner_id,
        kilobytes_used => $self->kilobytes_used,
    );
    $params{mount_path} = $self->mount_path if defined $self->mount_path;

    my $allocation = Genome::Disk::Allocation->create(%params);
    unless ($allocation) {
        Carp::confess "Could not create allocation!";
    }

    $self->status_message("Created allocation " . $allocation->id . " at " . $allocation->absolute_path);
    return 1;
}

1;

