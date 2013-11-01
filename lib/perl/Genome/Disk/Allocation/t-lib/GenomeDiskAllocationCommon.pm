package GenomeDiskAllocationCommon;
use base 'Test::Builder::Module';

use Exporter 'import';
our @EXPORT_OK = qw(create_test_volumes);

use Genome;

use File::Temp qw(tempdir);
use Filesys::Df qw();

sub create_test_volumes {
    my $n_pairs = shift;

    my $tb = __PACKAGE__->builder;

    my @volumes;
    $tb->subtest('create_test_volumes', sub {
        my $active_prefix = tempdir(
            'allocation_testing_active_XXXXXX',
            TMPDIR => 1,
            UNLINK => 1,
            CLEANUP => 1,
        );
        my $archive_prefix = tempdir(
            'allocation_testing_archive_XXXXXX',
            TMPDIR => 1,
            UNLINK => 1,
            CLEANUP => 1,
        );
        no warnings 'redefine';
        *Genome::Disk::Volume::active_volume_prefix  = sub { return $active_prefix };
        *Genome::Disk::Volume::archive_volume_prefix = sub { return $archive_prefix };
        use warnings;

        $Genome::Disk::Allocation::CREATE_DUMMY_VOLUMES_FOR_TESTING = 0;
        no warnings 'redefine';
        *Genome::Sys::current_user_is_admin = sub { return 1 };
        use warnings;

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

        # Make two dummy volumes, one to create allocation on and one to archive it to
        for my $n (1..$n_pairs) {
            for my $volume_prefix ($active_prefix, $archive_prefix) {
                my $volume_path = File::Spec->join($volume_prefix, 'volume' . $n);

                mkdir $volume_path
                    or die "mkdir failed: $!";

                use Genome::Disk::Volume;
                my $volume = Genome::Disk::Volume->create(
                    hostname => 'foo',
                    physical_path => 'foo/bar',
                    total_kb => Filesys::Df::df($volume_path)->{blocks},
                    disk_status => 'active',
                    can_allocate => '1',
                    mount_path => $volume_path,
                );
                $tb->ok($volume, 'made testing volume') or die;

                my $disk_group_name = (
                    $volume_prefix eq $active_prefix
                    ? 'testing_group'
                    : 'archive'
                );
                my $group = Genome::Disk::Group->get(disk_group_name => $disk_group_name);
                my $assignment = Genome::Disk::Assignment->create(
                    dg_id => $group->id,
                    dv_id => $volume->id,
                );
                $tb->ok($assignment, "assigned volume to group ($disk_group_name)");
                system("mkdir -p " . join('/', $volume->mount_path, $volume->groups->subdirectory));

                push @volumes, $volume;
            }

            $tb->is_num($volumes[-2]->is_archive, 0, 'first volume is not archive');
            $tb->is_num($volumes[-1]->is_archive, 1, 'second volume is archive');
        }
    });

    return @volumes;
}

1;
