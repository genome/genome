use strict;
use warnings;

use Test::More;

use Time::HiRes qw(usleep);
use above 'Genome';
use Genome::Disk::Allocation::Test qw(create_group create_tmpfs_volume create_barrier spawn_child waitpids);

# soft_limit_kb_tight_fit.t
# The purpose of this test is ...

my $group   = create_group('info_apipe');
for my $total_kb (600, 600, 1200) {
    my $volume = create_tmpfs_volume(total_kb => $total_kb, group => $group);
    ok($volume, sprintf('created %s', $volume->mount_path));
}
my $kilobytes_requested = 500;
my $max_n = 5;
my $barrier = create_barrier();
my @pids;
diag 'Errors below about not creating allocations are to be expected.';
for my $n (1..$max_n) {
    my $allocation_path = sprintf('test/soft_limit_kb_tight_fit/%s', $n); # this has to be three deep?
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

my @allocations = Genome::Disk::Allocation->get(disk_group_name => $group->disk_group_name);
is(scalar(@allocations), 4, 'got 4 allocations');

my @volumes = Genome::Disk::Volume->get(disk_group_names => $group->disk_group_name);
is(scalar(@volumes), 3, 'got 3 volumes');
for my $volume (@volumes) {
    ok($volume->allocated_kb <= $volume->total_kb,      sprintf('%s: allocated_kb <= total_kb'     , $volume->mount_path));
    ok($volume->allocated_kb <= $volume->hard_limit_kb, sprintf('%s: allocated_kb <= hard_limit_kb', $volume->mount_path));
    ok($volume->allocated_kb <= $volume->soft_limit_kb, sprintf('%s: allocated_kb <= soft_limit_kb', $volume->mount_path));
}

done_testing();
