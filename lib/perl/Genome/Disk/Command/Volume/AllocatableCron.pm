use strict;
use warnings;

package Genome::Disk::Command::Volume::AllocatableCron;

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
    return grep {   $_->is_near_soft_limit } @volumes;
}

sub select_volumes_to_enable {
    my $self = shift;
    my @volumes = @_;
    return grep { ! $_->is_near_soft_limit } @volumes;
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
        my $m = sprintf('%s: enabling', $v->mount_path);
        $self->$message_method($m);
        if (!$self->quiet) {
            if (! $v->is_used_over_soft_limit) {
                $self->$message_method("\tVolume is not used over soft limit.");
            }
            if (! $v->is_allocated_over_soft_limit) {
                $self->$message_method("\tVolume is not allocated over soft limit.");
            }
        }
        unless ($self->dry_run) {
            $v->can_allocate(1);
        }
    }
}
