package Genome::Disk::Command::Volume::Create;

use strict;
use warnings;

use Genome;

class Genome::Disk::Command::Volume::Create {
    is => 'Command::V2',
    doc => 'Create a new disk volume within a group',
    has_input => [
        mount_path => {
            is => 'DirectoryPath',
            doc => 'location where the volume is mounted on the filesystem',
        },
        total_kb => {
            is => 'Number',
            doc => 'Total size of the disk in kilobytes (or total size to make allocatable)',
        },
        hostname => {
            is => 'Text',
            doc => 'host of the storage volume being mounted',
        },
        physical_path => {
            is => 'Text',
            doc => 'path to the physical mount (in case the mount path is a logical volume thereunder)'
        },
        disk_group => {
            is => 'Genome::Disk::Group',
            doc => 'the group to which to initially assign the new volume',
        },
    ],
};

sub _is_hidden_in_docs { !Genome::Sys->current_user_is_admin }

sub help_brief {
    return 'create a new disk volume';
}

sub help_synopsis {
    return 'genome disk volume create --mount-path /path/to/mount --total-kb 1024 --hostname fileserver.example.com --phsyical-path /vol/to/mount --disk-group some_existing_group';
}

sub help_detail {
    return 'Creates a new disk volume within a group.  The volume needs to be currently mounted.';
}

sub execute {
    my $self = shift;

    my $mount_path = $self->mount_path;
    Genome::Sys->validate_existing_directory($mount_path);

    my $total_kb = int($self->total_kb);
    unless ($total_kb > 0) {
        $self->fatal_message('total_kb must be a positive integer');
    }

    my $hostname = $self->hostname;
    my $physical_path = $self->physical_path;
    my $disk_group = $self->disk_group;

    my $new_volume = Genome::Disk::Volume->create(
        mount_path => $mount_path,
        total_kb => $total_kb,
        hostname => $hostname,
        physical_path => $physical_path,
        disk_status => 'active',
        can_allocate => 1,
    );

    my $assignment = Genome::Disk::Assignment->create(
        volume_id => $new_volume->id,
        group_id => $disk_group->id,
    );

    return $new_volume->id;
}

sub exit_code_for_return_value {
    my $self = shift;
    my $return_value = shift;
    if (! $return_value) {
        $return_value = 1;
    } else {
        $return_value = 0;
    }
    return $return_value;
}

1;
