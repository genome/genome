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
    doc => 'move allocations from one volume to another',
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
    return 'moves allocations from one volume to another';
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    my $target_group = $self->target_group;
    my $target_volume = $self->target_volume;

    # If group and volume are specified, make sure the volume is in the group.
    return unless $target_group and $target_volume;
    return if $target_volume->groups(disk_group_name => $target_group->disk_group_name);
    push @errors, UR::Object::Tag->create(
        type => 'error',
        properties => [qw/ target_group target_volume /],
        desc => sprintf(
            'Specified target volume (%s) does not belong to group: %s',
            $target_volume->id, $target_group->disk_group_name
        ), 
    );

    return @errors;
}

sub execute {
    my $self = shift;

    for my $allocation ($self->allocations) {
        $allocation->move( $self->_resolve_move_params_for_allocation($allocation) );
    }

    $self->status_message("Successfully moved allocations!");
    return 1;
}

sub _resolve_move_params_for_allocation {
    my ($self, $allocation) = @_;

    my %params;
    if ($self->target_volume) {
        $params{target_mount_path} = $self->target_volume->mount_path;
    }

    if ($self->target_group) {
        $params{disk_group_name} = $self->target_group->disk_group_name;
    }

    if ( not $params{target_mount_path} and not $params{disk_group_name} ) {
        # Use the existing allocation's disk group
        $params{disk_group_name} = $allocation->disk_group_name;
    }

    return %params
}

1;

