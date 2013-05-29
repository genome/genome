package Genome::Disk::Command::Allocation::Move;

use strict;
use warnings;
use Genome;

class Genome::Disk::Command::Allocation::Move {
    is => 'Command::V2',
    has => [
        allocations => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'Allocations to move',
        }
    ],
    has_optional => [
        target_volume => {
            is => 'Genome::Disk::Volume',
            doc => 'Volume that allocations are to be moved to',
        },
        target_group => {
            is => 'Genome::Disk::Group',
            doc => 'Group that allocations are to be moved to',
        },
    ],
    doc => 'move alloations from one volume to another',
};

sub help_detail {
    return <<EOS
Moves allocations from one volume to another. Can specify a specific volume
to move them to or provide a group from which a volume will be selected.

If no volume or group is specified, then it will use the allocation's
disk_group_name as the target_group, the result being that the allocation
will be moved to another volume in the same group with enough free space to
hold it.
EOS
}

sub help_brief {
    return 'moves alloations from one volume to another';
}

sub execute {
    my $self = shift;

    for my $allocation ($self->allocations) {
        my %params;
        if ($self->target_volume) {
            $params{target_mount_path} = $self->target_volume->mount_path;
        }

        if ($self->target_group) {
            $params{disk_group_name} = $self->target_group->disk_group_name;
        } else {
            $self->status_message("Allocation " . $allocation->id . " has target_group ".$allocation->disk_group_name);
            $params{disk_group_name} = $allocation->disk_group_name;
        }

        $allocation->move(%params);
    }

    $self->status_message("Successfully moved allocations!");
    return 1;
}

1;

