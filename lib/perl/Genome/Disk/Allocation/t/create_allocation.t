use strict;
use warnings;

use Test::More;

use above 'Genome';
use Genome::Disk::Allocation::Test qw(create_group create_tmpfs_volume create_barrier spawn_child waitpids);

my $group_name = 'info_apipe';
my $group   = create_group($group_name);
my $volume  = create_tmpfs_volume(total_kb => 500, group => $group);
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
    my $allocation = Genome::Disk::Allocation->create(
        disk_group_name => $group->disk_group_name,
        allocation_path => $allocation_path,
        kilobytes_requested => 100,
        owner_class_name => 'UR::Value::Text',
        owner_id => 'Pwner',
    );
    is($allocation->allocation_path, $allocation_path, "create allocation ($allocation_path)");
}
done_testing();
