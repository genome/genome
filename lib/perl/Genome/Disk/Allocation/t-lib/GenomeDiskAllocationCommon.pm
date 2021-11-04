package GenomeDiskAllocationCommon;
use base 'Test::Builder::Module';

use strict;
use warnings;

use Exporter 'import';
our @EXPORT_OK = qw(create_test_volumes);

use Genome;

use File::Temp qw(tempdir);
use Filesys::Df qw();

sub create_test_volumes {
    my $n_volumes = shift;

    my $original_use_dummy_ids = UR::DataSource->use_dummy_autogenerated_ids;
    my $guard = Scope::Guard->new(sub { UR::DataSource->use_dummy_autogenerated_ids($original_use_dummy_ids) });
    UR::DataSource->use_dummy_autogenerated_ids(1);

    my $tb = __PACKAGE__->builder;

    my @volumes;
    $tb->subtest('create_test_volumes', sub {
        my $active_prefix = tempdir(
            'allocation_testing_active_XXXXXX',
            TMPDIR => 1,
            UNLINK => 1,
            CLEANUP => 1,
        );

        $Genome::Disk::Allocation::CREATE_DUMMY_VOLUMES_FOR_TESTING = 0;
        no warnings 'redefine';
        *Genome::Sys::current_user_is_admin = sub { return 1 };
        use warnings;

        # Make dummy group
        my $disk_group_name = 'testing_group';
        my $group = Genome::Disk::Group->create(
            disk_group_name => $disk_group_name,
            permissions => '755',
            setgid => '1',
            subdirectory => 'testing',
            unix_uid => 0,
            unix_gid => 0,
        );
        $tb->ok($group, 'successfully made testing group') or die;

        # Make a dummy volume for the allocation
        for my $n (1..$n_volumes) {
                my $volume_path = File::Spec->join($active_prefix, 'volume' . $n, 'Active');
                my $archive_path = File::Spec->join($active_prefix, 'volume' . $n, 'Archive');

                Genome::Sys->create_directory($volume_path);
                Genome::Sys->create_directory($archive_path);

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

                my $group = Genome::Disk::Group->get(disk_group_name => $disk_group_name);
                my $assignment = Genome::Disk::Assignment->create(
                    group_id => $group->id,
                    volume_id => $volume->id,
                );
                $tb->ok($assignment, "assigned volume to group ($disk_group_name)");
                system("mkdir -p " . join('/', $volume->mount_path, $volume->groups->subdirectory));

                push @volumes, $volume;
        }
    });

    return @volumes;
}

1;
