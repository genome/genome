package Genome::Disk::Volume;

use strict;
use warnings;

use Genome;
use Carp;

use Data::Dumper;
use Filesys::Df qw();
use List::Util qw(max);

class Genome::Disk::Volume {
    table_name => 'disk.volume',
    id_by => [
        dv_id => {is => 'Number', column_name => 'id' },
    ],
    has => [
        hostname => { is => 'Text' },
        physical_path => { is => 'Text' },
        mount_path => { is => 'Text' },
        disk_status => {
            is => 'Text',
            valid_values => ['inactive', 'active'],
        },
        can_allocate => {
            is => 'Number',
            valid_values => [0, 1],
        },

        # TODO remove this field when we switch to postgres
        _placeholder_creation_event_id => {
            is => 'Number',
            is_optional => 1,
            default_value => '1',
            column_name => 'creation_event_id',
            doc => 'This is here because of RT #88076.  Remove when we switch to postgres',
        },

        total_kb => { is => 'Number' },
        total_gb => {
            calculate_from => 'total_kb',
            calculate => q{ return int($total_kb / (2**20)) },
        },
        soft_limit_kb => {
            calculate_from => ['total_kb', 'maximum_reserve_size'],
            calculate => q{ $self->_compute_lower_limit($total_kb, 0.95, $maximum_reserve_size); },
        },
        soft_limit_gb => {
            calculate_from => 'soft_limit_kb',
            calculate => q{ return int($soft_limit_kb / (2**20)) },
        },
        hard_limit_kb => {
            calculate_from => ['total_kb', 'maximum_reserve_size'],
            calculate => q{ $self->_compute_lower_limit($total_kb, 0.98, int(0.5 * $maximum_reserve_size)); },
        },
        hard_limit_gb => {
            calculate_from => 'hard_limit_kb',
            calculate => q{ return int($hard_limit_kb / (2**20)) },
        },
        unallocated_kb => {
            calculate_from => ['total_kb', 'allocated_kb'],
            calculate => q{ return $total_kb - $allocated_kb },
        },
        cached_unallocated_kb => {
            is => 'Number',
            default_value => 0,
            column_name => 'unallocated_kb',
        },
        unallocated_gb => {
            calculate_from => 'unallocated_kb',
            calculate => q{ return int($unallocated_kb / (2**20)) },
        },
        allocated_kb => {
            is_calculated => 1,
        },
        percent_allocated => {
            calculate_from => ['total_kb', 'allocated_kb'],
            calculate => q{ return unless defined($total_kb) and defined($allocated_kb); return sprintf("%.2f", ( $allocated_kb / $total_kb ) * 100); },
        },
        used_kb => {
            calculate_from => ['mount_path'],
            calculate => q{ $self = shift; return $self->df->{used} },
        },
        percent_used => {
            calculate_from => ['total_kb', 'used_kb'],
            calculate => q{ return unless defined($used_kb); return sprintf("%.2f", ( $used_kb / $total_kb ) * 100); },
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
    id_generator => '-uuid',
    doc => 'Represents a particular disk volume (eg, sata483)',
};

sub mount {
    my $self = shift;
    return 0 if $self->is_mounted;
    my $mount_path = $self->mount_path;
    my $rel_path = $mount_path;
    $rel_path =~ s|/opt/gms/GMS1/||;
    unless (-d $mount_path) {
      Genome::Sys->create_directory($mount_path);
    }
    my $cmd = "curlftpfs ftp://clinus234/$rel_path $mount_path";
    Genome::Sys->shellcmd(cmd => $cmd);
    return 1;
}

sub unmount {
    my $self = shift;
    my $mount_path = $self->mount_path;
    return 0 if $self->is_local_volume;
    return 0 unless $self->is_mounted;
    my $cmd = "fusermount -u '$mount_path'";
    Genome::Sys->shellcmd(cmd => $cmd);
    return 1;
}

sub __display_name__ {
    my $self = shift;
    return $self->mount_path;
}

sub _compute_lower_limit {
    my $class = shift;
    my ($total_kb, $fraction, $maximum_reserve_size) = @_;
    my $fractional_limit = int($total_kb * $fraction);
    my $subtractive_limit = $total_kb - $maximum_reserve_size;
    return max($fractional_limit, $subtractive_limit);
}

sub allocated_kb {
    my $self = shift;

    # This used to use a copy of UR::Object::Set's server-side aggregate
    # logic to ensure no client-side calculation was done (for performance).
    # Now we're going to just let the set do the sum and monitor the
    # performance (in allocated_kb calculated property). If performance
    # is bad we can revert.

    my $set = Genome::Disk::Allocation->define_set(mount_path => $self->mount_path);
    my $field = 'kilobytes_requested';
    my $f = "sum($field)";

    # UR caches the value so we're just going to reach in and "fix" it.
    # Newer UR Sets store all their cached aggregate values under the __aggregates key
    if (exists $set->{__aggregates}) {
        $set->__invalidate_cache__($f);
    } elsif (exists $set->{$f}) {
        delete $set->{$f}
    }

    # We only want to time the sum aggregate, not including any time to connect to the DB
    Genome::Disk::Allocation->__meta__->data_source->get_default_handle();
    my $allocated_kb;
    Genome::Utility::Instrumentation::timer(
        'disk.volume.allocated_kb',
        sub { $allocated_kb = ($set->sum($field) or 0); }
    );

    # Now we'll check that it is cached so we test that the underlying
    # structure hasn't changed.
    unless(exists($set->{$f}) || (exists($set->{__aggregates}) && exists($set->{__aggregates}->{$f}))) {
        die $self->error_message("$f value not found in set's hash. Did underlying object structure change?");
    }

    return $allocated_kb;
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

our @dummy_volumes;
sub create_dummy_volume {
    my ($class, %params) = @_;
    my $mount_path = $params{mount_path};
    my $volume;
    if (!$mount_path || ($mount_path && $mount_path !~ /^\/tmp\//)) {
        $params{mount_path} = File::Temp::tempdir( 'tempXXXXX', TMPDIR => 1, CLEANUP => 1 );
        $volume = Genome::Disk::Volume->__define__(
            mount_path => $params{mount_path},
            total_kb => Filesys::Df::df($params{mount_path})->{blocks},,
            can_allocate => 1,
            disk_status => 'active',
            hostname => 'localhost',
            physical_path => '/tmp',
        );
        push @dummy_volumes, $volume;
        my $disk_group = Genome::Disk::Group->get(disk_group_name => $params{disk_group_name});
        Genome::Disk::Assignment->__define__(
            volume => $volume,
            group => $disk_group,
        );
    }
    else {
        $volume = Genome::Disk::Volume->get_active_volume(mount_path => $mount_path);
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

    if ($ENV{UR_DBI_NO_COMMIT}) {
        return 1;
    }

    # We can't use Filesys::Df::df because it doesn't report the mount path only the stats.
    my $mount_path = $self->mount_path;
    my @df_output = qx(df -P $mount_path 2> /dev/null);
    if ($! && $! !~ /No such file or directory/) {
        die $self->error_message(sprintf('Failed to `df %s` to check if volume is mounted: %s', $mount_path, $!));
    }

    return 1 if $self->is_local_volume;

    # any other case must be mounted somewhere deeper under / 
    my ($df_output) = grep { /\s$mount_path$/ } @df_output;
    return ($df_output ? 1 : 0);
}

sub is_local_volume {
    my $self = shift;
    my $mount_path = $self->mount_path;

    # consider the local disk of system first node to be mounted if mounted at /
    my $home = $ENV{GENOME_HOME};
    my $sysid = $ENV{GENOME_SYS_ID};
    my $dir = "$home/fs/$sysid";
    return ($mount_path =~ m|^$home| ? 1 : 0);
}
    


sub df {
    my $self = shift;

    unless ($self->is_mounted) {
        die $self->error_message(sprintf('Volume %s is not mounted!', $self->mount_path));
    }

    return Filesys::Df::df($self->mount_path);
}

sub sync_args {
    my $self = shift;
    my %args = @_;
    my $verbose = delete $args{verbose};
    if (keys %args) {
        die $self->error_message('Unexpected args to sync_meta: ' . join(', ', keys %args));
    }
    return (verbose => $verbose);
}

sub sync_usage {
    my $self = shift;
    my %args = $self->sync_args(@_);
    $self->sync_total_kb(%args);
    $self->sync_unallocated_kb(%args);
    return 1;
}

sub sync_total_kb {
    my $self = shift;
    my %args = $self->sync_args(@_);

    my $total_kb = $self->df->{blocks};
    if ($self->total_kb != $total_kb) {
        if ($args{verbose}) { $self->status_message(sprintf('Changing total_kb from %d to %d.', $self->total_kb, $total_kb)) }
        $self->total_kb($total_kb);
    }
}

sub sync_unallocated_kb {
    my $self = shift;
    my %args = $self->sync_args(@_);

    my $unallocated_kb = $self->unallocated_kb;
    if ($self->cached_unallocated_kb != $unallocated_kb) {
        if ($args{verbose}) { $self->status_message(sprintf('Changing cached_unallocated_kb from %d to %d.', $self->cached_unallocated_kb, $unallocated_kb)) }
        $self->cached_unallocated_kb($unallocated_kb);
    }
}

sub is_allocated_over_soft_limit {
    my $self = shift;
    return ($self->allocated_kb > $self->soft_limit_kb);
}

sub is_used_over_soft_limit {
    my $self = shift;
    return ($self->used_kb > $self->soft_limit_kb);
}

sub is_over_soft_limit {
    my $self = shift;
    return ($self->is_allocated_over_soft_limit || $self->is_used_over_soft_limit);
}

sub is_allocated_over_hard_limit {
    my $self = shift;
    return ($self->allocated_kb > $self->hard_limit_kb);
}

sub is_used_over_hard_limit {
    my $self = shift;
    return ($self->used_kb > $self->hard_limit_kb);
}

sub is_over_hard_limit {
    my $self = shift;
    return ($self->is_allocated_over_hard_limit || $self->is_used_over_hard_limit);
}

sub get_active_volume {
    my $class = shift;
    my %defaults = (disk_status => 'active', can_allocate => 1);
    my %params = (%defaults, @_);
    return $class->get(%params);
}

sub has_space {
    my ($self, $kilobytes_requested) = @_;

    return ($self->allocated_kb + $kilobytes_requested <= $self->soft_limit_kb);
}

1;

