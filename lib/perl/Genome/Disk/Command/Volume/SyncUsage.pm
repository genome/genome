package Genome::Disk::Command::Volume::SyncUsage;

use strict;
use warnings;

use Genome;

class Genome::Disk::Command::Volume::SyncUsage {
    is => 'Command::V2',
    has => {
        volumes => {
            is => 'Genome::Disk::Volume',
            is_many => 1,
            shell_args_position => 1,
            doc => 'Filter expression for volume(s) to sync.',
        },
        total_kb => {
            is => 'Boolean',
            default => 1,
            doc => 'Sync total_kb?',
        },
        unallocated_kb => {
            is => 'Boolean',
            default => 1,
            doc => 'Sync unallocated_kb?',
        },
    },
    doc => 'Sync total and unallocated KB usage for volumes',
};

sub help_detail {
    return __PACKAGE__->__meta__->doc;
}

sub execute {
    my $self = shift;
    $self->status_message('Disk volume sync usage...');

    for my $volume ( $self->volumes ) {
        $self->status_message('Syncing volume: %s', $volume->mount_path);
        for my $type (qw/ total unallocated /) {
            my $method = $type.'_kb';
            next unless $self->$method;
            my $kb = $volume->$method;
            my $sync_method = 'sync_'.$method;
            my $new_kb = $volume->$sync_method;
            $self->status_message('Synced %s kB from %d to %d', $type, $kb, $new_kb);
        }
    }

    $self->status_message('Disk volume sync usage...OK');
    return 1;
}

1;
