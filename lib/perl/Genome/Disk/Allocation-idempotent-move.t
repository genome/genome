#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 2;

use Genome::Disk::Allocation;

use File::Basename qw(dirname);

use lib File::Spec->join(dirname(__FILE__), 'Allocation', 't-lib');
use GenomeDiskAllocationCommon qw(create_test_volumes);

my @volumes = create_test_volumes(2);

subtest 'idempotent move' => sub {
    plan tests => 2;

    my $allocation_path = 'test_allocation_path';
    my $shadow_allocation_path =
        Genome::Disk::Detail::Allocation::Mover::_get_move_shadow_path(
            $allocation_path);


    my $allocation = Genome::Disk::Allocation->create(
        owner_class_name    => 'UR::Value',
        owner_id            => __FILE__,
        kilobytes_requested => 8,
        disk_group_name => $volumes[0]->disk_group_names,
        allocation_path => $allocation_path,
        mount_path => $volumes[0]->mount_path,
    );

    my $mover = Genome::Disk::Detail::Allocation::Mover->create(
        allocation_id => $allocation->id);
    my $shadow = Genome::Disk::Allocation->create(
        $mover->_get_move_shadow_params($allocation),
        disk_group_name => $volumes[2]->disk_group_names,
        mount_path => $volumes[2]->mount_path,
    );

    $allocation->move(
        disk_group_name => $volumes[2]->disk_group_names,
        target_mount_path => $volumes[2]->mount_path,
    );

    isnt($allocation->mount_path, $volumes[0]->mount_path, 'is not on old volume');
    is($allocation->mount_path, $volumes[2]->mount_path, 'is on new volume');
};
