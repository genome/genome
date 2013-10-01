package GenomeDiskAllocationCommon;
use base 'Test::Builder::Module';

use Exporter 'import';
our @EXPORT_OK = qw(create_test_volumes);

use Genome;

use File::Temp qw(tempdir);
use Filesys::Df qw();

sub create_test_volumes {
    my $tb = __PACKAGE__->builder;

    my $test_dir = tempdir(
        'allocation_testing_XXXXXX',
        TMPDIR => 1,
        UNLINK => 1,
        CLEANUP => 1,
    );

    $Genome::Disk::Allocation::CREATE_DUMMY_VOLUMES_FOR_TESTING = 0;
    no warnings 'redefine';
    *Genome::Sys::current_user_is_admin = sub { return 1 };
    use warnings;

    # Make two dummy volumes, one to create allocation on and one to archive it to
    my @volumes;
    for (1..2) {
        my $volume_path = tempdir(
            "test_volume_" . $_ . "_XXXXXXX",
            DIR => $test_dir,
            CLEANUP => 1,
            UNLINK => 1,
        );
        my $volume = Genome::Disk::Volume->create(
            id => $_,
            hostname => 'foo',
            physical_path => 'foo/bar',
            mount_path => $volume_path,
            total_kb => Filesys::Df::df($volume_path)->{blocks},
            disk_status => 'active',
            can_allocate => '1',
        );
        $tb->ok($volume, 'made testing volume') or die;
        push @volumes, $volume;
    }

    # Make dummy groups
    my $group = Genome::Disk::Group->create(
        disk_group_name => 'testing_group',
        permissions => '755',
        setgid => '1',
        subdirectory => 'testing',
        unix_uid => 0,
        unix_gid => 0,
    );
    $tb->ok($group, 'successfully made testing group') or die;

    my $archive_group = Genome::Disk::Group->create(
        disk_group_name => 'archive',
        permissions => '755',
        setgid => '1',
        subdirectory => 'testing',
        unix_uid => 0,
        unix_gid => 0,
    );
    $tb->ok($archive_group, 'created archive group');
    push @Genome::Disk::Allocation::APIPE_DISK_GROUPS, $group->disk_group_name, $archive_group->disk_group_name;

    # Assign volumes to groups
    my $assignment = Genome::Disk::Assignment->create(
        dg_id => $group->id,
        dv_id => $volumes[0]->id,
    );
    $tb->ok($assignment, 'assigned one volume to testing group');
    system("mkdir -p " . join('/', $volumes[0]->mount_path, $volumes[0]->groups->subdirectory));

    my $archive_assignment = Genome::Disk::Assignment->create(
        dg_id => $archive_group->id,
        dv_id => $volumes[1]->id,
    );
    $tb->ok($archive_assignment, 'created one volume in archive group');
    system("mkdir -p " . join('/', $volumes[1]->mount_path, $volumes[1]->groups->subdirectory));

    # Redefine archive and active volume prefix, make sure this make is_archive return expected values
    no warnings 'redefine';
    *Genome::Disk::Volume::active_volume_prefix = sub { return $volumes[0]->mount_path };
    *Genome::Disk::Volume::archive_volume_prefix = sub { return $volumes[1]->mount_path };
    use warnings;

    $tb->is_num($volumes[0]->is_archive, 0, 'first volume is not archive');
    $tb->is_num($volumes[1]->is_archive, 1, 'second volume is archive');
    $tb->is_eq($volumes[0]->archive_mount_path, $volumes[1]->mount_path, 'first volume correctly identifies its archive volume');

    return @volumes;
}

1;
