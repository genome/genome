package Genome::Disk::Command::Volume::AllocatableCron;

use strict;
use warnings;

use Genome;

use File::stat qw(stat);
use Genome::Utility::File::Mode qw(mode);

class Genome::Disk::Command::Volume::AllocatableCron {
    is => 'Command::V2',
    has => [
        dry_run => {
            is => 'Boolean',
            default => 0,
        },
        quiet => {
            is => 'Boolean',
            default => 0,
        },
        disable => {
            is => 'Boolean',
            default => 1,
        },
        enable => {
            is => 'Boolean',
            default => 1,
        },
        volumes => {
            is => 'Genome::Disk::Volume',
            is_many => 1,
            shell_args_position => 1,
            require_user_verify => 0,
        },
    ],
};

sub execute {
    my $self = shift;

    my @volumes = $self->volumes;
    my @active_volumes = grep { $_->disk_status eq 'active' } @volumes;

    if ($self->disable) {
        my @volumes_to_disable = $self->select_volumes_to_disable(@active_volumes);
        $self->disable_volumes(@volumes_to_disable);
    }

    if ($self->enable) {
        my @volumes_to_enable  = $self->select_volumes_to_enable(@active_volumes);
        $self->enable_volumes(@volumes_to_enable);
    }

    return 1;
}

sub print_message {
    my $self = shift;
    my $message = join('', @_);
    chomp $message;
    print $message, "\n";
}

sub select_volumes_to_disable {
    my $self = shift;
    my @volumes = @_;
    return grep { $_->can_allocate && $_->is_near_soft_limit } @volumes;
}

sub select_volumes_to_enable {
    my $self = shift;
    my @volumes = @_;
    return grep { ! $_->can_allocate && ! $_->is_near_soft_limit } @volumes;
}

sub disable_volumes {
    my $self = shift;
    my @volumes = @_;

    my $message_method = (($self->dry_run || !$self->quiet) ? 'print_message' : 'debug_message');

    for my $v (@volumes) {
        my $m = sprintf('%s: disabling', $v->mount_path);
        $self->$message_method($m);
        if (!$self->quiet) {
            if ($v->is_used_over_soft_limit) {
                $self->$message_method("\tVolume is used over soft limit.");
            }
            if ($v->is_allocated_over_soft_limit) {
                $self->$message_method("\tVolume is allocated over soft limit.");
            }
        }
        unless ($self->dry_run) {
            $v->can_allocate(0);
        }
    }
}

sub enable_volumes {
    my $self = shift;
    my @volumes = @_;

    my $message_method = (($self->dry_run || !$self->quiet) ? 'print_message' : 'debug_message');

    for my $v (@volumes) {
        my $error = sanity_check_for_enable($v);
        if ($error) {
            my $m = sprintf('%s: skipping: %s', $v->mount_path, $error);
            $self->$message_method($m);
            next;
        }
        else {
            my $m = sprintf('%s: enabling', $v->mount_path);
            $self->$message_method($m);
        }

        unless ($self->dry_run) {
            $v->can_allocate(1);
        }
    }
}

sub sanity_check_for_enable {
    my $volume = shift;

    unless (-d $volume->mount_path) {
        return 'mount_path is not a directory';
    }

    my $disk_group = $volume->groups;

    my $subdir_path = File::Spec->join($volume->mount_path, $disk_group->subdirectory);
    unless (-d $subdir_path) {
        return 'expected subdirectory is not a directory';
    }

    my $subdir_stat = stat($subdir_path);
    unless ($subdir_stat->gid == $disk_group->unix_gid) {
        return sprintf 'subdirectory is not owned by %s', $disk_group->group_name;
    }

    my $subdir_mode = mode($subdir_path);
    unless ($subdir_mode->is_setgid) {
        return 'subdirectory is not setgid';
    }
    unless ($subdir_mode->is_group_rwx) {
        return 'subdirectory is not rwx';
    }

    if ($volume->is_used_over_soft_limit) {
        return 'volume usage is over soft limit';
    }

    if ($volume->is_allocated_over_soft_limit) {
        return 'volume is allocated over soft limit';
    }

    return;
}
