#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;
use Test::MockObject;

plan tests => 3;

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

subtest 'succesful sync' => sub{
    plan tests => 3;

    my $cmd = $class->execute(
        volumes => [$volume],
    );
    ok($cmd->result, 'execute');
    $volume->called_ok('sync_total_kb');
    $volume->called_ok('sync_unallocated_kb');
    $volume->clear;

};

my $unmounted_volume = Test::MockObject->new();
$unmounted_volume->set_isa('Genome::Disk::Volume');
$unmounted_volume->set_always('mount_path', '/unmounted/001');
$unmounted_volume->set_always('total_kb', 20);
my $err = 'Volume /unmounted_volume/001 is not mounted!';
$unmounted_volume->mock('sync_total_kb', sub{ die $err });
$unmounted_volume->set_always('unallocated_kb', 300);
$unmounted_volume->set_always('sync_unallocated_kb', 300);

subtest 'successful despite unmounted volume' => sub{
    plan tests => 6;

    my $cmd = $class->execute(
        volumes => [$volume, $unmounted_volume],
    );
    ok($cmd->result, 'execute');
    like($cmd->error_message, qr/$err/, 'correct error message');
    $volume->called_ok('sync_total_kb');
    $volume->called_ok('sync_unallocated_kb');
    $volume->clear;
    $unmounted_volume->called_ok('sync_total_kb');
    $unmounted_volume->called_ok('sync_unallocated_kb');
    $unmounted_volume->clear;

};

done_testing();
