package Genome::Disk::Volume;

use strict;
use warnings;

use Genome;
use Carp;

class Genome::Disk::Volume {
    table_name => 'DISK_VOLUME',
    id_by => [
        dv_id => { is => 'Number', column_name => 'id'},
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
        unallocated_kb => { is => 'Number' },
        unallocated_gb => {
            calculate_from => 'unallocated_kb',
            calculate => q{ return int($unallocated_kb / (2**20)) },
        },
        allocated_kb => { 
            calculate_from => ['total_kb','unallocated_kb'],
            calculate => q{ return ($total_kb - $unallocated_kb); },
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
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'Represents a particular disk volume (eg, sata483)',
};

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
        $params{mount_path} = File::Temp::tempdir( TEMPLATE => 'tempXXXXX', CLEANUP => 1 );
        $volume = Genome::Disk::Volume->__define__(
            mount_path => $params{mount_path},
            unallocated_kb => 104857600, # 100 GB
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

1;

