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

use above 'Genome';
use GenomeDiskAllocationTest qw(create_tmpfs_volume create_barrier spawn_child waitpids);

my $volume  = create_tmpfs_volume(total_kb => 500);
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
    my ($group) = $volume->groups;
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
