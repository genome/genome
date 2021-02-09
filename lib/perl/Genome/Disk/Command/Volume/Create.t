#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Test::More tests => 10;
use Test::Exception;

my $class = 'Genome::Disk::Command::Volume::Create';

use_ok($class);

my $mount_path = Genome::Sys->create_temp_directory();
my $total_kb = 1024;
my $hostname = Genome::Sys->hostname();
my $physical_path = $mount_path;

my $disk_group = Genome::Disk::Group->__define__(
    name => 'disk_command_volume_t_test_group' . time(),
    id => (0 - time()),
);

my $bad_mount_path = '/dev/null/bonus/path';

my $cmd1 = $class->create(
    mount_path => $bad_mount_path,
    total_kb => $total_kb,
    hostname => $hostname,
    physical_path => $physical_path,
    disk_group => $disk_group,
);
isa_ok($cmd1, $class, 'created first command');
throws_ok(sub { $cmd1->execute }, qr/directory/, 'fails with bad directory');

my $cmd2 = $class->create(
    mount_path => $mount_path,
    total_kb => 0,
    hostname => $hostname,
    physical_path => $physical_path,
    disk_group => $disk_group,
);
isa_ok($cmd2, $class, 'created second command');
throws_ok(sub { $cmd2->execute }, qr/must be a positive integer/, 'fails with bad total_kb');

my $cmd3 = $class->create(
    mount_path => $mount_path,
    total_kb => $total_kb,
    hostname => $hostname,
    physical_path => $physical_path,
    disk_group => $disk_group,
);
isa_ok($cmd3, $class, 'created third command');
my $volume_id = $cmd3->execute;
ok($volume_id, 'command executes');

my $volume = Genome::Disk::Volume->get(mount_path => $mount_path);
isa_ok($volume, 'Genome::Disk::Volume', 'got new volume');
is($volume->id, $volume_id, 'command returned new volume id');

is($volume->disk_group_names, $disk_group->name, 'volume is assigned to group');

