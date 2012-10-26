package Genome::Disk::Volume;

use strict;
use warnings;

use Genome;
use Carp;

use Data::Dumper;
use Filesys::Df qw();

class Genome::Disk::Volume {
    table_name => 'DISK_VOLUME',
    id_by => [
        dv_id => {is => 'Number'},
    ],
    has => [
        hostname => { is => 'Text' },
        physical_path => { is => 'Text' },
        mount_path => { is => 'Text' },
        disk_status => { is => 'Text' },
        can_allocate => { is => 'Number' },

        total_kb => { is => 'Number' },
        total_gb => {
            calculate_from => 'total_kb',
            calculate => q{ return int($total_kb / (2**20)) },
        },

        soft_limit_kb => {
            calculate_from => 'total_kb',
            calculate => q{ return int($total_kb * 0.95); },
        },
        soft_limit_gb => {
            calculate_from => 'soft_limit_kb',
            calculate => q{ return int($soft_limit_kb / (2**20)) },
        },
        hard_limit_kb => {
            calculate_from => 'total_kb',
            calculate => q{ return int($total_kb * 0.98); },
        },
        hard_limit_gb => {
            calculate_from => 'hard_limit_kb',
            calculate => q{ return int($hard_limit_kb / (2**20)) },
        },

        unallocated_kb => {
            calculate_from => ['total_kb', 'allocated_kb'],
            calculate => q{ $self = shift;
                            Carp::cluck('Cannot set calculated unallocated_kb') if (@_);
                            return $total_kb - $allocated_kb;
            },
        },
        old_unallocated_kb => {
            is => 'Number',
            default_value => 0,
            column_name => 'unallocated_kb',
        },
        unallocated_gb => {
            calculate_from => 'unallocated_kb',
            calculate => q{ return int($unallocated_kb / (2**20)) },
        },

        allocated_kb => {
            calculate_from => ['mount_path'],
            calculate => q/
                my $meta        = Genome::Disk::Allocation->__meta__;
                my $table_name  = $meta->table_name;
                my $data_source = $meta->data_source;
                my $owner       = $data_source->owner;

                my $query_string;
                if ($data_source->isa('UR::DataSource::Oracle')) {
                    my $fq_table_name = join('.', $owner, $table_name);
                    $query_string = sprintf(q(select sum(kilobytes_requested) from %s where mount_path = ?), $fq_table_name);
                } elsif ($data_source->isa('UR::DataSource::Pg') || $data_source->isa('UR::DataSource::SQLite')) {
                    $query_string = sprintf(q(select sum(kilobytes_requested) from %s where mount_path = ?), $table_name);
                } else {
                    die sprintf('allocated_kb cannot be calculated for %s', $data_source->class);
                }
                my $dbh = $data_source->get_default_handle();
                my $query_handle = $dbh->prepare($query_string);
                $query_handle->bind_param(1, $mount_path);
                $query_handle->execute();

                my $row_arrayref = $query_handle->fetchrow_arrayref();
                $query_handle->finish();
                unless (defined $row_arrayref) {
                    Genome::Disk::Volume->error_message("Could not calculate allocated kb from database.");
                }
                unless (1 == scalar(@$row_arrayref)) {
                    Genome::Disk::Volume->error_message(sprintf(
                        "Incorrect number of elements returned from SQL query (%s) expected 1.",
                        scalar(@$row_arrayref))
                    );
                }

                return ($row_arrayref->[0] or 0);
            /,
        },
        old_allocated_kb => {
            calculate_from => ['total_kb', 'old_unallocated_kb'],
            calculate => q{ return $total_kb - $old_unallocated_kb; },
        },
        percent_allocated => {
            calculate_from => ['total_kb', 'allocated_kb'],
            calculate => q{ return sprintf("%.2f", ( $allocated_kb / $total_kb ) * 100); },
        },

        used_kb => {
            calculate_from => ['mount_path'],
            calculate => sub {
                my $mount_path = shift;
                return 0 unless -e $mount_path;
                my ($used_kb) = qx(df -Pk $mount_path | grep $mount_path | awk '{print \$3}') =~ /(\d+)/; 
                return $used_kb
            },
        },
        percent_used => {
            calculate_from => ['total_kb', 'used_kb'],
            calculate => q{ return sprintf("%.2f", ( $used_kb / $total_kb ) * 100); },
        },

        unallocatable_reserve_size => {
            calculate_from => ['total_kb', 'unallocatable_volume_percent', 'maximum_reserve_size'],
            calculate => q{
                my $buffer = int($total_kb * $unallocatable_volume_percent);
                $buffer = $maximum_reserve_size if $buffer > $maximum_reserve_size;
                return $buffer;
            },
            doc => 'Size of reserve in kb that cannot be allocated to but can still be used by reallocations',
        },
        unusable_reserve_size => {
            calculate_from => ['total_kb', 'unusable_volume_percent', 'unallocatable_reserve_size'],
            calculate => q{
                my $buffer = int($total_kb * $unusable_volume_percent);
                $buffer = $unallocatable_reserve_size if $buffer > $unallocatable_reserve_size;
                return $buffer;
            },
            doc => 'Size of reserve in kb that cannot be allocated to in any way',
        },
        allocatable_kb => {
            calculate_from => ['unallocated_kb', 'unallocatable_reserve_size'],
            calculate => q{ 
                my $allocatable = $unallocated_kb - $unallocatable_reserve_size;
                $allocatable = 0 if $allocatable < 0;  # Possible due to reallocation having a smaller reserve
                return $allocatable;
            },
        },
        allocatable_gb => {
            calculate_from => 'allocatable_kb',
            calculate => q{ return int($allocatable_kb / (2**20)) },
        },

        unallocatable_volume_percent => {
            is => 'Number',
            is_constant => 1,
            is_classwide => 1,
            column_name => '',
            value => '.05',
        },
        unusable_volume_percent => {
            is => 'Number',
            is_constant => 1,
            is_classwide => 1,
            column_name => '',
            value => '.02',
        },

        maximum_reserve_size => {
            is => 'Number',
            is_constant => 1,
            is_classwide => 1,
            column_name => '',
            value => 1_073_741_824,
        },
    ],

    has_many_optional => [
        disk_group_names => {
            via => 'groups',
            to => 'disk_group_name',
        },
        groups => {
            is => 'Genome::Disk::Group',
            via => 'assignments',
            to =>  'group',
        },
        assignments => {
            is => 'Genome::Disk::Assignment',
            reverse_id_by => 'volume',
        },
        allocations => {
            is => 'Genome::Disk::Allocation',
            calculate_from => 'mount_path',
            calculate => q| return Genome::Disk::Allocation->get(mount_path => $mount_path); |,
        },
    ],
    data_source => 'Genome::DataSource::Oltp',
    doc => 'Represents a particular disk volume (eg, sata483)',
};

sub __display_name__ {
    my $self = shift;
    return $self->mount_path;
}

sub get_lock {
    my ($class, $mount_path, $tries) = @_;
    $tries ||= 120;
    my $modified_mount = $mount_path;
    $modified_mount =~ s/\//_/g;
    my $volume_lock = Genome::Sys->lock_resource(
        resource_lock => $ENV{GENOME_LOCK_DIR} . '/allocation/volume' . $modified_mount,
        max_try => $tries,
        block_sleep => 1,
        wait_announce_interval => 10,
    );
    return $volume_lock;
}

sub create_dummy_volume {
    my ($class, %params) = @_;
    my $mount_path = $params{mount_path};
    my $volume;
    if (!$mount_path || ($mount_path && $mount_path !~ /^\/tmp\//)) {
        $params{mount_path} = File::Temp::tempdir( 'tempXXXXX', TMPDIR => 1, CLEANUP => 1 );
        $volume = Genome::Disk::Volume->__define__(
            mount_path => $params{mount_path},
            total_kb => 104857600,
            can_allocate => 1,
            disk_status => 'active',
            hostname => 'localhost',
            physical_path => '/tmp',
        );
        my $disk_group = Genome::Disk::Group->get(disk_group_name => $params{disk_group_name});
        Genome::Disk::Assignment->__define__(
            volume => $volume,
            group => $disk_group,
        );
    }
    else {
        $volume = Genome::Disk::Volume->get(mount_path => $mount_path, disk_status => 'active', can_allocate => 1);
    }

    return $volume;
}

# FIXME Rather than having a mix of archive/active volume logic here, it would probably be cleaner 
# to have archive/active subclasses of volume. Would need a subclassify method, but that's it.
sub archive_volume_prefix {
    return '/gscarchive';
}

sub active_volume_prefix {
    return '/gscmnt';
}

sub is_archive {
    my $self = shift;
    my $archive_prefix = $self->archive_volume_prefix;
    if ($self->mount_path =~ /^$archive_prefix/) {
        return 1;
    }
    return 0;
}

sub archive_mount_path {
    my $self = shift;
    return if $self->is_archive;
    my $mount = $self->mount_path;
    my $archive_volume_prefix = $self->archive_volume_prefix;
    my $active_volume_prefix = $self->active_volume_prefix;
    $mount =~ s/$active_volume_prefix/$archive_volume_prefix/;
    return $mount;
}

sub archive_volume {
    my $self = shift;
    return if $self->is_archive;
    return Genome::Disk::Volume->get(mount_path => $self->archive_mount_path);
}

sub active_mount_path {
    my $self = shift;
    return if not $self->is_archive;
    my $mount = $self->mount_path;
    my $archive_volume_prefix = $self->archive_volume_prefix;
    my $active_volume_prefix = $self->active_volume_prefix;
    $mount =~ s/$archive_volume_prefix/$active_volume_prefix/;
    return $mount;
}

sub active_volume {
    my $self = shift;
    return if not $self->is_archive;
    return Genome::Disk::Volume->get(mount_path => $self->active_mount_path);
}

sub is_mounted {
    my $self = shift;

    # We can't use Filesys::Df::df because it doesn't report the mount path only the stats.
    my $mount_path = $self->mount_path;
    my @df_output = qx(df -P $mount_path 2> /dev/null);
    if ($! && $! !~ /No such file or directory/) {
        die $self->error_message(sprintf('Failed to `df %s` to check if volume is mounted: %s', $mount_path, $!));
    }

    my ($df_output) = grep { /\s$mount_path$/ } @df_output;
    return ($df_output ? 1 : 0);
}

sub df {
    my $self = shift;

    unless ($self->is_mounted) {
        die $self->error_message(sprintf('Volume %s is not mounted!', $self->mount_path));
    }

    return Filesys::Df::df($self->mount_path);
}

sub sync_usage {
    my $self = shift;
    my %args = @_;
    my $verbose = delete $args{verbose};
    if (keys %args) {
        die $self->error_message('Unexpected args to sync_meta: ' . join(', ', keys %args));
    }

    my $total_kb = $self->df->{blocks};
    my $unallocated_kb = $self->unallocated_kb;

    if ($self->total_kb != $total_kb) {
        if ($verbose) { $self->status_message(sprintf('Changing total_kb from %d to %d.', $self->total_kb, $total_kb)) }
        $self->total_kb($total_kb);
    }
    if ($self->old_unallocated_kb != $unallocated_kb) {
        if ($verbose) { $self->status_message(sprintf('Changing old_unallocated_kb from %d to %d.', $self->old_unallocated_kb, $unallocated_kb)) }
        $self->old_unallocated_kb($unallocated_kb);
    }

    return 1;
}

sub is_over_soft_limit {
    my $self = shift;
    return (
        $self->allocated_kb > $self->soft_limit_kb
        || $self->used_kb > $self->soft_limit_kb
    );
}

sub is_over_hard_limit {
    my $self = shift;
    return (
        $self->allocated_kb > $self->hard_limit_kb
        || $self->used_kb > $self->hard_limit_kb
    );
}

1;

