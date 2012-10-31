use strict;
use warnings;

use Test::More;

use above 'Genome';
use Genome::Disk::Allocation::Test qw(create_tmpfs_volume create_barrier spawn_child waitpids);

my $volume  = create_tmpfs_volume(total_kb => 500);
my ($group) = $volume->groups;

my $file_kb = $volume->soft_limit_kb - $volume->used_kb + 1;
my $filename = join('/', $volume->mount_path, 'space_filler');
system("dd if=/dev/zero of=${filename} bs=1k count=${file_kb}") && die "dd failed: $!";
is(-s $filename, $file_kb * 1024, "created ${file_kb} kB file") or die;
diag "used_kb = " . $volume->used_kb . "\n";
diag "soft_limit_kb = " . $volume->soft_limit_kb . "\n";
ok($volume->is_over_soft_limit, 'volume should now be over soft limit') or die;

my $allocation_path = 'test/create_allocation/some_allocation'; # this has to be three deep?
{
    my @child_allocations = Genome::Disk::Allocation->get('allocation_path like' => "$allocation_path/%");
    is(scalar(@child_allocations), 0, 'no child allocations');
}
{
    my $allocation = Genome::Disk::Allocation->get(allocation_path => $allocation_path);
    is($allocation, undef, 'allocation does not already exist');
}
{
    diag 'Errors below about not creating allocations are to be expected.';
    my $allocation = eval { Genome::Disk::Allocation->create(
        disk_group_name => $group->disk_group_name,
        allocation_path => $allocation_path,
        kilobytes_requested => 100,
        owner_class_name => 'UR::Value::Text',
        owner_id => 'Pwner',
    )};
    diag 'Errors above about not creating allocations are to be expected.';
    is($allocation, undef, 'no allocation should have been created');
}
done_testing();
