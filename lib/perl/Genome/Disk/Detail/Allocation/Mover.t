#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::Exception;
use Test::More;

use_ok('Genome::Disk::Detail::Allocation::Mover') or die;

subtest 'mover fails when target group and volume are specified' => sub {
    plan tests => 1;

    throws_ok(
        sub{
            Genome::Disk::Detail::Allocation::Mover->create(
                allocation_id=> 'allocation',
                disk_group_name => 'group',
                target_mount_path => 'volume'
            );
        },
        qr/Can not specify both target group and volume/, 
        'move fails when suppling target group and volume',
    );

};

done_testing();
