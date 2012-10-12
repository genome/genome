package Genome::Disk::Allocation::Test;
@ISA = qw(Exporter);
@EXPORT_OK = qw(create_group create_tmpfs_volume create_barrier spawn_child waitpids);

use strict;
use warnings;

use Carp qw(croak);
use Test::More;
use Time::HiRes qw(usleep);
use above 'Genome';

BEGIN {
    if ( ! -t STDIN ) {
        plan skip_all => 'This test can only be run interactively!';
    }
    if ( $ENV{UR_DBI_NO_COMMIT} ) {
        plan skip_all => 'This test can not be run with UR_DBI_NO_COMMIT enabled!';
    }
};

our $imported = 0;
sub import {
    unless ($imported) {
        protect_real_data_sources();
        if ($ENV{UR_DBI_NO_COMMIT}) {
            die 'This is mean to be run with UR_DBI_NO_COMMIT not enabled.';
        }

        my $lds = Genome::DataSource::LocalDataSource->get();
        $lds->hijack_class('Genome::Disk::Assignment') or die;
        $lds->hijack_class('Genome::Disk::Group') or die;
        $lds->hijack_class('Genome::Disk::Volume') or die;
        $lds->hijack_class('Genome::Disk::Allocation') or die;

    $Genome::Disk::Allocation::_execute_system_command_template =
        '%s %s -e "use above Genome; use Genome::Disk::Allocation::Test; %s->%s(%s); UR::Context->commit;"';
        $imported = 1;
    }

    # Now call Exporter's import, a.k.a. export_to_level...
    Genome::Disk::Allocation::Test->export_to_level(1, @_);
}

sub create_tmpfs {
    our %tmpfs_mnt_pts;
    my %arg = @_;
    my $size = delete $arg{size} || die;
    my $path = delete $arg{path};

    unless ($path) {
        $path = File::Temp::tempdir('tmpfs_XXXX', TMPDIR => 1);
    }

    # TODO: check exit code
    system("sudo mount -t tmpfs -o size=$size tmpfs $path");

    push @{$tmpfs_mnt_pts{$$}}, $path;

    return $path;
}

sub create_tmpfs_volume {
    my %arg = @_;
    my $total_kb = delete $arg{total_kb} || die;
    my $group = delete $arg{group} || die;

    my $mount_path = create_tmpfs(size => "${total_kb}k");

    my $volume = Genome::Disk::Volume->create(
        disk_status    => 'active',
        hostname       => 'localhost',
        physical_path  => $mount_path,
        mount_path     => $mount_path,
        total_kb       => $total_kb,
        unallocated_kb => $total_kb,
        can_allocate   => 1,
    );

    my $assignment = Genome::Disk::Assignment->create(
        dg_id => $group->id,
        dv_id => $volume->id,
    );

    UR::Context->commit;

    return $volume;
}

sub protect_real_data_sources {
    my @ds_module_filenames = grep { /Genome\/DataSource/ } keys %INC;
    if (@ds_module_filenames) {
        warn "Genome::DataSources already loaded:\n" . join("\n", @ds_module_filenames);
        die('Genome::DataSources already loaded! This test runs without UR_DBI_NO_COMMIT enabled so we do not allow this.');
    }

    my $lds = Genome::DataSource::LocalDataSource->get();

    my %ds_module_filename = (
        'Genome::DataSource::Main' => 'Genome/DataSource/Main.pm',
        'Genome::DataSource::GMSchema' => 'Genome/DataSource/GMSchema.pm',
    );
    for my $ds_module (keys %ds_module_filename) {
        my $ds_module_filename = $ds_module_filename{$ds_module};
        $INC{$ds_module_filename} = $lds->__meta__->module_path;
        UR::Object::Type->define(
            class_name => $ds_module,
            is => $lds->class,
        );
    }
}

sub create_group {
    my $disk_group_name = shift;
    my $group = Genome::Disk::Group->create(
        disk_group_name => $disk_group_name,
        permissions => '755',
        sticky => '1',
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
    for my $path (@{$tmpfs_mnt_pts{$$}}) {
        system("sudo umount $path");
        if (-d $path) {
            system("rmdir $path");
        }
    }
};

1;
