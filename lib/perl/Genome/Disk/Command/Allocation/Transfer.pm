package Genome::Disk::Command::Allocation::Transfer;

use strict;
use warnings;

use Genome;
use Genome::Disk::Allocation;

use File::Spec;

class Genome::Disk::Command::Allocation::Transfer {
    is => 'Command::V2',
    has => [
        allocations => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'allocations to transfer',
        },
        target_group => {
            is => 'Genome::Disk::Group',
            doc => 'remote group to which to transfer the allocation',
        },
        target_volume => {
            is => 'Genome::Disk::Volume',
            doc => 'remote volume to which to transfer the allocation',
        },
        remote_host => {
            is => 'Text',
            doc => 'ssh-enabled remote host for transfer',
        },
        remote_user => {
            is => 'Text',
            doc => 'remote username for ssh',
        },
        remote_port => {
            is => 'Number',
            doc => 'remote port to use for SSH',
            default_value => 22,
        },
        authorization_key => {
            is => 'Text',
            doc => 'path to the SSH key to use for rsync transfer',
        },
    ],
};

sub help_detail {
    return <<EOHELP
Transfers an allocation to a remote volume.

This command requires an SSH server to be running on a host connected to the remote volume.  It will rsync the allocation to the remote volume and delete the contents from the local volume.
EOHELP

}

sub help_brief {
    return 'transfers allocations to a remote volume';
}

sub execute {
    my $self = shift;

    my $volume = $self->target_volume;
    my $group = $self->target_group;

    my $assignment = Genome::Disk::Assignment->get(group_id => $group->id, volume_id => $volume->id);
    unless ($assignment) {
        $self->fatal_message('Volume %s is not a member of Group %s', $volume->mount_path, $group->name);
    }

    for my $allocation ($self->allocations) {
        $self->status_message('Starting %s...', $allocation->__display_name__);

        #It's customary to shell out for allocation operations so they are committed immediately and independently.
        require Genome::Disk::Allocation;
        Genome::Disk::Allocation::_execute_system_command(
            $self->class,
            '_process_allocation',
            id => $allocation->id,
            group_id => $group->id,
            volume_id => $volume->id,
            host => $self->remote_host,
            user => $self->remote_user,
            key => $self->authorization_key,
            port => $self->remote_port,
        );

        $self->status_message('Finished %s.', $allocation->__display_name__);
    }

    return 1;
}

sub _process_allocation {
    my $class = shift;
    my %params = @_;

    my ($allocation_object, $allocation_lock) = Genome::Disk::Allocation->get_with_lock($params{id});
    my $target_volume = Genome::Disk::Volume->get($params{volume_id});
    my $target_group = Genome::Disk::Group->get($params{group_id});

    my $host = $params{host};
    my $port = $params{port};
    my $user = $params{user};
    my $key  = $params{key};

    $class->transfer($allocation_object, $target_volume, $target_group, $host, $port, $user, $key);

    #clean up the old directory on commit
    $allocation_object->_create_observer(
        $allocation_object->_get_deletion_observers
    );

    $class->update_allocation($allocation_object, $target_volume, $target_group);

    Genome::Disk::Allocation->_commit_unless_testing();
    $allocation_lock->unlock();
}

sub _reload_allocation {
    my $class = shift;

    Genome::Disk::Allocation->_reload_allocation(@_);
}

sub transfer {
    my $class = shift;
    my $allocation = shift;
    my $target_volume = shift;
    my $target_group = shift;
    my $host = shift;
    my $port = shift;
    my $user = shift;
    my $key = shift;

    my $source_path = $allocation->absolute_path .'/';

    my $destination_path = File::Spec->join(
        $target_volume->mount_path,
        $target_group->subdirectory,
        $allocation->allocation_path,
    );

    my @ssh_opts = ('ssh', '-p', $port, '-i', $key, '-o', 'StrictHostKeyChecking=no');

    my $remote_destination = sprintf('%s@%s:%s', $user, $host, $destination_path);

    #first ensure directory exists as target for rsync
    Genome::Sys->shellcmd(
        cmd => [@ssh_opts, sprintf('%s@%s', $user, $host), 'mkdir', '-p', $destination_path],
    );
    Genome::Sys->shellcmd(
        cmd => ['rsync', '-av', '-e', join(' ', @ssh_opts), $source_path, $remote_destination],
    );
}

sub update_allocation {
    my $class = shift;
    my $allocation = shift;
    my $target_volume = shift;
    my $target_group = shift;

    $allocation->mount_path($target_volume->mount_path);
    $allocation->group_subdirectory($target_group->subdirectory);
    $allocation->disk_group_name($target_group->name);

    $allocation->_update_owner_for_move;
}

1;
