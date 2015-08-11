package Genome::Disk::Command::Volume::SyncUsage;

use strict;
use warnings;

use Genome;

use Try::Tiny;

class Genome::Disk::Command::Volume::SyncUsage {
    is => ['Genome::Role::Logger', 'Command::V2'],
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
    has_transient => {
        tie_stderr => {
            is => 'Boolean',
            default => 1,
            doc => '(warning) globally tie STDERR to this logger',
        },
    },
    doc => 'Sync total and unallocated KB usage for volumes',
};

sub help_detail {
    return __PACKAGE__->__meta__->doc;
}

sub execute {
    my $self = shift;

    for my $volume ( $self->volumes ) {
        my $transaction = UR::Context::Transaction->begin;
        try {
            $self->info('Syncing volume: '.$volume->mount_path);
            if ($self->total_kb)      {
                $self->debug('Sync total kB...');
                $volume->sync_total_kb;
            }
            if ($self->unallocated_kb) {
                $self->debug('Sync unallocated kB...');
                $volume->sync_unallocated_kb;
            }
            $transaction->commit or die 'Failed to commit volume! '.$volume->mount_path;
        }
        catch {
            $self->error($_);
            $transaction->rollback;
        };
    }

    return 1;
}

1;
