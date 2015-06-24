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
    doc => 'Sync usage info for volume (e.g. total KB and unallocated KB)',
};

sub help_detail {
    'Sync usage info for volume (e.g. total KB and unallocated KB).'
}

sub execute {
    my $self = shift;

    my %args = (verbose => 1);
    for my $volume ( $self->volumes ) {
        my $transaction = UR::Context::Transaction->begin;
        try {
            $self->info('Syncing volume: '.$volume->mount_path);
            if ($self->total_kb)      {
                $self->debug('Sync total kB...');
                $volume->sync_total_kb(%args);
            }
            if ($self->unallocated_kb) {
                $self->debug('Sync unallocated kB...');
                $volume->sync_unallocated_kb(%args);
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
