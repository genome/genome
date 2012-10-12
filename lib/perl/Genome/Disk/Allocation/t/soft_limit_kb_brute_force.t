use strict;
use warnings;

use Test::More;

use Time::HiRes qw(usleep);
use above 'Genome';
use Genome::Disk::Allocation::Test qw(create_group create_tmpfs_volume create_barrier spawn_child waitpids);

# soft_limit_kb.t
# The purpose of this test is to determine whether the soft_limit_kb is respected. To do that we attempt to create many more allocations than a volume can hold and see how full the volume gets.

my $total_kb = 1000;
my $kilobytes_requested = 50;
my $max_n = int($total_kb / $kilobytes_requested) + 5;

my $group_name = 'info_apipe';
my $group   = create_group($group_name);
my $volume  = create_tmpfs_volume(total_kb => $total_kb, group => $group);
my $barrier = create_barrier();
my @pids;
diag 'Errors below about not creating allocations are to be expected.';
for my $n (1..$max_n) {
    my $allocation_path = sprintf('testing/%s', $n);
    push @pids, spawn_child(
        barrier => $barrier,
        closure => sub {
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
        },
    );
}
usleep(500_000);
unlink $barrier;
waitpids(@pids);
diag 'Errors above about not creating allocations are to be expected.';

ok($volume->allocated_kb <= $volume->total_kb,      sprintf('%s: allocated_kb <= total_kb'     , $volume->mount_path));
ok($volume->allocated_kb <= $volume->hard_limit_kb, sprintf('%s: allocated_kb <= hard_limit_kb', $volume->mount_path));
ok($volume->allocated_kb <= $volume->soft_limit_kb, sprintf('%s: allocated_kb <= soft_limit_kb', $volume->mount_path));

my $tolerance_kb = $volume->soft_limit_kb - $kilobytes_requested;
ok($volume->allocated_kb > $tolerance_kb, sprintf('%s: allocated_kb > %d', $volume->mount_path, $tolerance_kb))
    or diag "allocated_kb = " . $volume->allocated_kb;

done_testing();
