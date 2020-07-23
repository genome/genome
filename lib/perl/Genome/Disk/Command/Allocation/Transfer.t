#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Genome::Test::Factory::DiskAllocation;
use Sub::Install qw(reinstall_sub);
use Test::More tests => 7;

my $class = 'Genome::Disk::Command::Allocation::Transfer';
use_ok($class);

my $owner = UR::Value::Text->get('dummy-owner');
my $allocation = Genome::Test::Factory::DiskAllocation->setup_object(
    owner => $owner,
);

my $destination_group = Genome::Disk::Group->create(
    disk_group_name => 'testing_group_for_transfer_t',
    permissions => '755',
    setgid => '1',
    subdirectory => 'testing123',
    unix_uid => 0,
    unix_gid => 0,
);

my $vol_path = Genome::Sys->create_temp_directory;
my $destination_volume = Genome::Disk::Volume->create(
    hostname => 'foo',
    physical_path => '/bar',
    mount_path => $vol_path,
    total_kb => 100_000,
    disk_status => 'active',
    can_allocate => 1,
);

Genome::Disk::Assignment->create(
    volume_id => $destination_volume->id,
    group_id => $destination_group->id,
);


my $cmd = $class->create(
    allocations => [$allocation],
    target_group => $destination_group,
    target_volume => $destination_volume,
    remote_host => 'localhost', #this won't be a real test since we have no remote server
    remote_user => Genome::Sys->current_user->username,
    #remote_port => 22, #default
    authorization_key => 'dummy',
);
isa_ok($cmd, $class, 'instantiated class');

reinstall_sub({
    into => 'Genome::Sys',
    as => 'shellcmd',
    code => sub {
        $class->debug_message('Skipping shellcmd execution for testing');
        return 1;
    },
});

ok($cmd->execute, 'executed command');

is($allocation->disk_group_name, $destination_group->disk_group_name, 'group updated');
is($allocation->mount_path, $destination_volume->mount_path, 'volume updated');
is($allocation->group_subdirectory, $destination_group->subdirectory, 'subdirectory updated');

SKIP: {
    skip 'shellcmd would be required to rsync data', 1;
    ok(-e $allocation->absolute_path, 'new location was created');
};

done_testing();
