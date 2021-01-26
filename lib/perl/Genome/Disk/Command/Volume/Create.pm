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

1;
