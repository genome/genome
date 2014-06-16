#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above "Genome";

require File::Temp;
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
is($obj->disk_allocation, undef, 'no disk_allocations, yet');
is($obj->disk_allocations, undef, 'no disk_allocations, yet');

my $tempdir = File::Temp::tempdir(CLEANUP => 1);
my @volumes = (
    Genome::Disk::Volume->__define__(
        mount_path => $tempdir,
        disk_status => 'active',
        can_allocate => 1,
        total_kb => 10000,
    ),
    Genome::Disk::Volume->__define__(
        mount_path => Genome::Disk::Volume->archive_volume_prefix,
        disk_status => 'active',
        can_allocate => 1,
        total_kb => 10000,
    ),
);
ok(!$volumes[0]->is_archive, 'volume 1 is not archived');
ok($volumes[1]->is_archive, 'volume 2 is archived');

my @expected_allocations = Genome::Disk::Allocation->__define__(
    owner_id => $obj->id,
    owner_class_name => $obj->class,
    mount_path => $volumes[0]->mount_path,
    disk_group_name => 'info_apipe',
    group_subdirectory => '',
    allocation_path => '1',
    kilobytes_requested => '2000',
);
ok(@expected_allocations, 'create allocation 1');
mkdir $expected_allocations[0]->absolute_path;

my $disk_allocation = $obj->disk_allocations;
is_deeply($disk_allocation, $expected_allocations[0], 'disk_allocation');
my @disk_allocations = $obj->disk_allocations;
is_deeply(\@disk_allocations, \@expected_allocations, 'disk_allocations');

push @expected_allocations, Genome::Disk::Allocation->__define__(
    owner_id => $obj->id,
    owner_class_name => $obj->class,
    mount_path => $volumes[0]->mount_path,
    disk_group_name => 'info_apipe',
    group_subdirectory => '',
    allocation_path => '1',
    kilobytes_requested => '2000',
);
is(@expected_allocations, 2, 'create allocation 2');
mkdir $expected_allocations[1]->absolute_path;

my $sorter = sub{ [ sort { $a cmp $b } @_ ]; }; # sort the mem locations :)
is(eval{ $obj->disk_allocation; }, undef, 'singular disk_allocation called when there are more than one');
@disk_allocations = $obj->disk_allocations;
is_deeply($sorter->(@disk_allocations), $sorter->(@expected_allocations), 'disk_allocations');

my @associated_disk_allocations = $obj->associated_disk_allocations;
is_deeply($sorter->(@associated_disk_allocations), $sorter->(@expected_allocations), 'associated_disk_allocations');

# reallocate
is(List::Util::sum(map { $_->kilobytes_requested } @expected_allocations), 4000, 'kilobytes_requested');
ok($obj->reallocate, 'reallocate');
isnt(List::Util::sum(map { $_->kilobytes_requested } @expected_allocations), 4000, 'kilobytes_requested after reallocate');

# archivable
is_deeply([ map { $_->archivable } @expected_allocations ], [qw/ 1 1 /], 'expected allocation are archivable');
is($obj->archivable, 1, 'archivable');
$expected_allocations[1]->archivable(0);
is_deeply([ map { $_->archivable } @expected_allocations ], [qw/ 1 0 /], 'expected allocation 2 is not archivable');
is($obj->archivable, '','not archivable');
$expected_allocations[1]->archivable(1);

# is_archived
is($obj->is_archived, '', 'not archived');
$expected_allocations[1]->mount_path( $volumes[1]->mount_path );
is($obj->is_archived, 1,'archived');

# unarchive
Sub::Install::reinstall_sub({
        code => sub{ $_[0]->mount_path( $volumes[0]->mount_path ); return 1; },
        into => 'Genome::Disk::Allocation',
        as   => 'unarchive',
    });
ok($obj->unarchive, 'unarchive');
is($obj->is_archived, '', 'not archived');

# deallocate
ok($obj->deallocate, 'deallocate');
isa_ok($expected_allocations[0], 'UR::DeletedRef');
isa_ok($expected_allocations[1], 'UR::DeletedRef');

# delete
ok($obj->delete, 'delete');

#print join("\n", map { $_->absolute_path } @expected_allocations); print Data::Dumper::Dumper(\@expected_allocations); <STDIN>;
done_testing();
