#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above "Genome";

require Genome::Test::Factory::DiskAllocation;
require List::Util;
require Sub::Install;
use Test::More;

use_ok('Genome::Utility::ObjectWithAllocations') or die;

my $meta = UR::Object::Type->define(
    class_name => 'GenomeTest::HasAllocations',
    is => [qw/ Genome::Utility::ObjectWithAllocations /],
    #has => { dummy_val => { is => 'Text' }, },
);
ok($meta, 'defined has allocations test class');

my $obj = GenomeTest::HasAllocations->create();
ok($obj, 'create');
is($obj->disk_allocation, undef, 'no disk_allocation, yet');
is($obj->disk_allocations, undef, 'no disk_allocations, yet');

my @expected_allocations = Genome::Test::Factory::DiskAllocation->generate_obj(owner => $obj);
ok(@expected_allocations, 'create allocation 1');

my $disk_allocation = $obj->disk_allocations;
is_deeply($disk_allocation, $expected_allocations[0], 'disk_allocation');
my @disk_allocations = $obj->disk_allocations;
is_deeply(\@disk_allocations, \@expected_allocations, 'disk_allocations');

push @expected_allocations, Genome::Test::Factory::DiskAllocation->generate_obj(owner => $obj);
is(@expected_allocations, 2, 'create allocation 2');
mkdir $expected_allocations[1]->absolute_path;

my $sorter = sub{ [ sort { $a cmp $b } @_ ]; }; # sort the mem locations :)
is(eval{ $obj->disk_allocation; }, undef, 'singular disk_allocation called when there are more than one');
@disk_allocations = $obj->disk_allocations;
is_deeply($sorter->(@disk_allocations), $sorter->(@expected_allocations), 'disk_allocations');

my @associated_disk_allocations = $obj->associated_disk_allocations;
is_deeply($sorter->(@associated_disk_allocations), $sorter->(@expected_allocations), 'associated_disk_allocations');

# reallocate_disk_allocations
is(List::Util::sum(map { $_->kilobytes_requested } @expected_allocations), 4000, 'kilobytes_requested');
ok($obj->reallocate_disk_allocations, 'reallocate_disk_allocations');
isnt(List::Util::sum(map { $_->kilobytes_requested } @expected_allocations), 4000, 'kilobytes_requested after reallocate_disk_allocations');

# are_disk_allocations_archivable
is_deeply([ map { $_->archivable } @expected_allocations ], [qw/ 1 1 /], 'expected allocation are archivable');
is($obj->are_disk_allocations_archivable, 1, 'disk_allocations are archivable');
$expected_allocations[1]->archivable(0);
is_deeply([ map { $_->archivable } @expected_allocations ], [qw/ 1 0 /], 'expected allocation 2 is not archivable');
is($obj->are_disk_allocations_archivable, '','disk_allocations are not archivable');
$expected_allocations[1]->archivable(1);

# are_disk_allocations_archived
is($obj->are_disk_allocations_archived, '', 'disk_allocations are not archived');
$expected_allocations[1]->archive;
is($obj->are_disk_allocations_archived, 1, 'disk_allocations are archived');

# unarchive
ok($obj->unarchive_disk_allocations, 'unarchive_disk_allocations');
is($obj->are_disk_allocations_archived, '', 'disk allocations are not archived');

# deallocate_disk_allocations
ok($obj->deallocate_disk_allocations, 'deallocate_disk_allocation');
isa_ok($expected_allocations[0], 'UR::DeletedRef');
isa_ok($expected_allocations[1], 'UR::DeletedRef');

# delete
push @expected_allocations, Genome::Test::Factory::DiskAllocation->generate_obj(owner => $obj);
is(@expected_allocations, 3, 'create new allocation to test delete');
ok($obj->delete, 'delete');
ok(UR::Context->commit, 'commit');
isa_ok($obj, 'UR::DeletedRef');
isa_ok($expected_allocations[2], 'UR::DeletedRef');

#print join("\n", map { $_->absolute_path } @expected_allocations); print Data::Dumper::Dumper(\@expected_allocations); <STDIN>;
done_testing();
