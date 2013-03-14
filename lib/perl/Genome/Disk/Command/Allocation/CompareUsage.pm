package Genome::Disk::Command::Allocation::CompareUsage;

use strict;
use warnings;

use Genome;

use Data::Dumper;
use IO::File qw();
use IO::Handle qw();
use Text::CSV qw();

class Genome::Disk::Command::Allocation::CompareUsage {
    is => 'Command::V2',
    has => [
        allocations => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'Allocations for which to compare disk usage with allocation size',
        },
    ],

    has_optional => [
        output => {
            is => 'Text',
            doc => 'Output filename, defaults to STDOUT',
        },

        du => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Use `du` to get actual disk usage (more accurate, but slower)',
        },

        update_cache => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Update kilobytes_used field on each allocation (implies --du)',
        },

        commit_chunk_size => {
            is => 'Integer',
            default_value => 1000,
            doc => 'Maximum number of cache updates to commit at once',
        },
    ],
    doc => 'Dumps a csv comparing `du` output with kilobytes_requested for each allocation',
};


sub help_synopsis {
    return "    Compare kilobytes_requested with kilobytes_used for one or more allocations.\n";
}

sub help_detail {
    return <<EOS;
Output format is CSV, with columns:

ID,CREATION_TIME,REALLOCATION_TIME,KILOBYTES_REQUESTED,KILOBYTES_USED
EOS
}


sub execute {
    my $self = shift;

    my $file = $self->_open_file();
    my $csv = $self->_get_csv($file);

    my $uncommited_count = 0;

    for my $a ($self->allocations) {
        my $stats = $self->_collect_stats($a);

        $self->_write_row($file, $csv, $stats);

        if ($self->update_cache) {
            $uncommited_count++;
            $self->_update_stats($a, $stats);

            if (($uncommited_count > $self->commit_chunk_size)) {
                UR::Context->commit();
            }
        }
    }

    if ($self->update_cache && $uncommited_count) {
        UR::Context->commit();
    }

    return scalar $self->allocations;
}

sub _open_file {
    my $self = shift;

    if ($self->output) {
        return IO::File->new($self->output, 'w');
    } else {
        return IO::Handle->new_from_fd(*STDOUT, 'w');
    }
}

sub _get_csv {
    my ($self, $file) = @_;

    my $csv = Text::CSV->new({binary => 1});

    return $csv;
}

sub _collect_stats {
    my ($self, $a) = @_;

    my $stats = {
        'id'                  => $a->id,
        'creation_time'       => $a->creation_time,
        'reallocation_time'   => $a->reallocation_time,
        'kilobytes_requested' => $a->kilobytes_requested,
        'kilobytes_used'      => $self->_get_kilobytes_used($a),
    };

    return $stats;
}

sub _get_kilobytes_used {
    my ($self, $a) = @_;

    if ($self->du) {
        return $a->du();
    } else {
        return $a->kilobytes_used;
    }
}

my @_DISPLAY_COLUMN_ORDER = (
    'id',
    'creation_time',
    'reallocation_time',
    'kilobytes_requested',
    'kilobytes_used',
);

sub _write_row {
    my ($self, $file, $csv, $stats) = @_;

    my @cols = map {$stats->{$_}} @_DISPLAY_COLUMN_ORDER;
    my $status = $csv->combine(@cols);
    # XXX check status for error
    $file->printf("%s\n", $csv->string);
}

sub _update_stats {
    my ($self, $a, $stats) = @_;

    $a->kilobytes_used($stats->{'kilobytes_used'});
}

1;
