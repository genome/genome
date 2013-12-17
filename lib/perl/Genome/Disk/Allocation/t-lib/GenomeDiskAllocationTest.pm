package GenomeDiskAllocationTest;
@ISA = qw(Exporter);
@EXPORT_OK = qw(create_group create_tmpfs_volume create_barrier spawn_child waitpids);

use strict;
use warnings;

use Carp qw(croak);
use File::Basename qw(dirname);
use Test::More skip_all => 'Test requires fuseext2 and leaves a mess that has to be cleaned up manually if correct permissions are not setup.';
use Time::HiRes qw(usleep);
use Sys::Hostname qw(hostname);
use above 'Genome';

BEGIN {
    if ( $ENV{UR_DBI_NO_COMMIT} ) {
        plan skip_all => 'This test can not be run with UR_DBI_NO_COMMIT enabled!';
    }
};

sub import {
    # Don't change @_!
    my $class = $_[0];

    protect_real_data_sources();

    require Genome::DataSource::LocalDataSource;
    Genome::DataSource::LocalDataSource->import(qw(
        Genome::Disk::Allocation
        Genome::Disk::Assignment
        Genome::Disk::Group
        Genome::Disk::Volume
    ));

    # ensure test environment is loaded for _execute_system_command
    my $lib_dir = dirname(__FILE__);
    for my $option ("-I$lib_dir", '-MGenomeDiskAllocationTest') {
        unless (grep { $_ eq $option } @Genome::Disk::Allocation::_execute_system_command_perl5opt) {
            push @Genome::Disk::Allocation::_execute_system_command_perl5opt, $option;
        }
    }

    $class->export_to_level(1, @_);
}


sub class_to_filename {
    my $class = shift;
    my $filename = $class;
    $filename =~ s/::/\//g;
    $filename .= '.pm';
    return $filename;
}

sub protect_real_data_sources {
    my $lds = Genome::DataSource::LocalDataSource->get();
    my $lds_path = $lds->__meta__->module_path;

    my %ds_module_filename = (
        'Genome::DataSource::GMSchema' => 'Genome/DataSource/GMSchema.pm',
    );
    for my $ds_module (keys %ds_module_filename) {
        my $ds_module_filename = $ds_module_filename{$ds_module};

        if ($INC{$ds_module_filename}) {
            if ($INC{$ds_module_filename} eq $lds_path) {
                warn "$ds_module already protected!";
            } else {
                die "$ds_module already loaded; cannot protect it!";
            }
        }

        $INC{$ds_module_filename} = $lds_path;
        UR::Object::Type->define(
            class_name => $ds_module,
            is => $lds->class,
        );
    }
}

sub create_tmpfs {
    our %tmpfs_mnt_pts;
    my %arg = @_;
    my $size = delete $arg{size} || die;
    my $mount_path = delete $arg{path};

    unless ($mount_path) {
        $mount_path = File::Temp::tempdir('allocation_test_mount_XXXX', TMPDIR => 1);
    }

    my (undef, $fs_path) = File::Temp::tempfile('allocation_test_fs_XXXX', TMPDIR => 1, OPEN => 1);
    system("dd if=/dev/zero of=${fs_path} bs=1k count=${size}") && die "dd failed: $!";
    system("/sbin/mkfs.ext2 -F ${fs_path}") && die "mkfs failed: $!";

    my $hostname = hostname();
    system("fuseext2 -o rw+ ${fs_path} ${mount_path}") && die sprintf('%s: fuseext2 failed: %s', $hostname, $!);

    push @{$tmpfs_mnt_pts{$$}}, [$fs_path, $mount_path];

    return $mount_path;
}

sub create_tmpfs_volume {
    my %arg = @_;
    my $total_kb = delete $arg{total_kb} || die;
    my $group = delete $arg{group} || Genome::Disk::Group->get(disk_group_name => $ENV{GENOME_DISK_GROUP_DEV});
    unless ($group) {
        $group = create_group($ENV{GENOME_DISK_GROUP_DEV});
    }

    my $mount_path = create_tmpfs(size => "${total_kb}");

    my $volume = Genome::Disk::Volume->create(
        disk_status    => 'active',
        hostname       => 'localhost',
        physical_path  => $mount_path,
        mount_path     => $mount_path,
        total_kb       => $total_kb,
        can_allocate   => 1,
    );

    my $assignment = Genome::Disk::Assignment->create(
        dg_id => $group->id,
        dv_id => $volume->id,
    );

    UR::Context->commit;

    return $volume;
}

sub create_group {
    my $disk_group_name = shift;
    my $group = Genome::Disk::Group->create(
        disk_group_name => $disk_group_name,
        permissions => '755',
        setgid => '1',
        subdirectory => $disk_group_name,
        unix_uid => 0,
        unix_gid => 0,
    );
    return $group;
}

sub create_barrier {
    my ($barrier_fh, $barrier) = File::Temp::tempfile(
        'genome_disk_allocation_barrier_XXXXXX',
        TMPDIR => 1,
        UNLINK => 0,
    );
    $barrier_fh->close();
    return $barrier;
}

sub spawn_child {
    my %args = @_;
    my $closure = delete $args{closure} || croak 'required argument missing: closure';
    my $usleep  = delete $args{usleep} || 10_000;
    my $barrier = delete $args{barrier};

    if (keys %args) {
        croak 'unexpected argument(s):', join(', ', keys %args);
    }

    my $pid = UR::Context::Process->fork();
    if (not defined $pid) {
        die 'Failed to fork()!';
    } elsif ($pid) {
        # Parent
        return $pid;
    } else {
        # Child
        while ($barrier && -e $barrier) {
            usleep($usleep);
        }
        $closure->();
        exit(0);
    }
}

sub waitpids {
    while (my $pid = shift @_) {
        my $wp_code = waitpid($pid, 0);
        if (-1 == $wp_code) {
            warn "Child pid already reaped.";
        } elsif (0 == $wp_code) {
            warn $?;
        }
    }
}

END {
    our %tmpfs_mnt_pts;
    for my $ref (@{$tmpfs_mnt_pts{$$}}) {
        my ($fs_path, $mount_path) = @$ref;

        system("fusermount -u ${mount_path}") && warn "fusermount failed";

        if (-d $mount_path) {
            system("find $mount_path -type d -delete") && warn "find -delete failed";
        }

        if (-f $fs_path) {
            system("rm ${fs_path}") && warn "rm failed";
        }
    }
};

1;
