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
$volume->set_true('sync_total_kb');
$volume->set_true('sync_unallocated_kb');

my $cmd = $class->execute(
    volumes => [$volume],
    tie_stderr => 0,
);
ok($cmd->result, 'execute');
$volume->called_ok('sync_total_kb');
$volume->called_ok('sync_unallocated_kb');

done_testing();
