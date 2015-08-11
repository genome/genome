#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;
use Test::MockObject;

my $class = 'Genome::Disk::Command::Volume::SyncUsage';
use_ok($class) or die;

my $volume = Test::MockObject->new();
$volume->set_always('mount_path', '/gsctmp/001');
$volume->set_true('is_mounted');
$volume->set_isa('Genome::Disk::Volume');
$volume->set_always('total_kb', 10);
$volume->set_always('sync_total_kb', 100);
$volume->set_always('unallocated_kb', 200);
$volume->set_always('sync_unallocated_kb', 200);

my $cmd = $class->execute(
    volumes => [$volume],
);
ok($cmd->result, 'execute');
$volume->called_ok('sync_total_kb');
$volume->called_ok('sync_unallocated_kb');

done_testing();
