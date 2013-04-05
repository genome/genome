use strict;
use warnings;

use Test::More;

use File::Basename qw(dirname);
our $lib_dir;
BEGIN {
    $lib_dir = dirname(__FILE__) . '-lib' ;
};
BEGIN {
    # Untested. Site::TGI changes broke these tests.
    use lib $lib_dir;
};

use Time::HiRes qw(usleep);
use above 'Genome';
use GenomeDiskAllocationTest qw(create_tmpfs_volume create_barrier spawn_child waitpids);

my $group_total_kb = 0;
my $group;
my $n_volumes = 10;
for my $n (1..$n_volumes) {
    my $total_kb = 200;
    $group_total_kb += $total_kb;
    my $volume = create_tmpfs_volume(total_kb => $total_kb);
    ($group) = $volume->groups;
    ok($volume, sprintf('created %s', $volume->mount_path));
}

my $allocated_kb = 0;
my $i = 1;
my $barrier = create_barrier();
my @pids;
diag 'Errors below about not creating allocations are to be expected.';
while ($allocated_kb < int(0.95 * $group_total_kb)) {
    my $allocations_per_child = 10;
    my $kilobytes_requested = 10 + int(rand(20));
    push @pids, spawn_child(
        barrier => $barrier,
        closure => sub {
            for my $j (1..$allocations_per_child) {
                my $allocation_path = sprintf('test/many_many/%d_%d', $i, $j);
                my $allocation = Genome::Disk::Allocation->create(
                    disk_group_name => $group->disk_group_name,
                    allocation_path => $allocation_path,
                    kilobytes_requested => $kilobytes_requested,
                    owner_class_name => 'UR::Value::Text',
                    owner_id => 'Pwner',
                );
                unless ($allocation) {
                    warn "Allocation ($allocation_path) not found!";
                }
            }
        },
    );
    $allocated_kb += $allocations_per_child * $kilobytes_requested;
    $i++;
}
usleep(500_000);
unlink $barrier;
waitpids(@pids);
diag 'Errors above about not creating allocations are to be expected.';

my @allocations = Genome::Disk::Allocation->get(disk_group_name => $group->disk_group_name);
ok(@allocations > 1, 'got some allocations');

my @volumes = Genome::Disk::Volume->get(disk_group_names => $group->disk_group_name);
is(scalar(@volumes), $n_volumes, "got $n_volumes volumes");
for my $volume (@volumes) {
    ok($volume->allocated_kb <= $volume->total_kb,            sprintf('%s: allocated_kb <= total_kb'     , $volume->mount_path));
    ok($volume->allocated_kb <= $volume->hard_limit_kb,       sprintf('%s: allocated_kb <= hard_limit_kb', $volume->mount_path));
    ok($volume->allocated_kb <= $volume->soft_limit_kb,       sprintf('%s: allocated_kb <= soft_limit_kb', $volume->mount_path));
    ok($volume->allocated_kb >= int(0.5 * $volume->total_kb), sprintf('%s: allocated_kb >= 50%% total_kb', $volume->mount_path));
}

done_testing();
