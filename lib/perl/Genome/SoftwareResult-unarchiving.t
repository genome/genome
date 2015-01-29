
use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above "Genome";
use Genome::Disk::Allocation;
use Genome::Test::Factory::SoftwareResult::User;

use Test::More tests => 8;

class Genome::Test_SR {
    is => 'Genome::SoftwareResult',
    has_param => [
        p => { is => 'Text' },
    ],
    has_input => [
        i => { is => 'Text' },
    ],
};

my $sr = Genome::Test_SR->create(p => 'test', i => 'test');
isa_ok($sr, 'Genome::SoftwareResult', 'setup a test software result');

my $disk_allocation = Genome::Disk::Allocation->create(owner_id => $sr->id, owner_class_name => $sr->class, kilobytes_requested => 1, disk_group_name => $ENV{GENOME_DISK_GROUP_MODELS}, allocation_path => $sr->id . '-alloc');

no warnings qw(redefine);
my $simulate_is_archived = 1;
sub Genome::Disk::Allocation::is_archived { return $simulate_is_archived; }

my $call_count = 0;
sub Genome::Disk::Allocation::unarchive { $call_count++; return 1; }
use warnings;

ok($disk_allocation->is_archived, 'setup an archived allocation for the test software result');

my $result_users = Genome::Test::Factory::SoftwareResult::User->setup_user_hash;

my $sr_again = Genome::Test_SR->get_with_lock(p => 'test', i => 'test', users => $result_users);
is($sr_again, $sr, 'got the same result again');
is($call_count, 1, 'the allocation was unarchived in get_with_lock');

my $sr_third = Genome::Test_SR->get_or_create(p => 'test', i => 'test', users => $result_users);
is($sr_third, $sr, 'got the same result a third time');
is($call_count, 2, 'the allocation was unarchived in get_or_create');

$simulate_is_archived = 0;
my $sr_fourth = Genome::Test_SR->get_with_lock(p => 'test', i => 'test', users => $result_users);
is($sr_fourth, $sr, 'got the same result a fourth time');
is($call_count, 2, 'the allocation was not unarchived in get_with_lock when not archived');

