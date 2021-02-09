#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Test::More tests => 9;
use Test::Exception;

my $class = 'Genome::Disk::Command::Group::Create';

use_ok($class);

my $test_name = 'disk_command_group_t_test_group' . time;

my $bad_subdir = 'dir1/dir2';
my $good_subdir = 'dir12';

my $unix_group = Genome::Config::get('sys_group');

my $cmd1 = $class->create(
    name => $test_name,
    subdirectory => $bad_subdir,
    unix_group => $unix_group,
);
isa_ok($cmd1, $class, 'created first command');
throws_ok(sub { $cmd1->execute }, qr/contains invalid characters/, 'fails with bad subdirectory');

my $cmd2 = $class->create(
    name => $test_name,
    subdirectory => $good_subdir,
    unix_group => $unix_group,
);
isa_ok($cmd2, $class, 'created second command');
my $disk_group_id = $cmd2->execute;
ok($disk_group_id, 'command executed successfully');

my $disk_group = Genome::Disk::Group->get(name => $test_name);
isa_ok($disk_group, 'Genome::Disk::Group', 'created disk group');
is($disk_group->id, $disk_group_id, 'command returned id of new group');

my $cmd3 = $class->create(
    name => $test_name,
    subdirectory => $good_subdir,
    unix_group => $unix_group,
);
isa_ok($cmd3, $class, 'created third command');
throws_ok(sub {$cmd3->execute }, qr/already exists/, 'fails to create new group with same name');


