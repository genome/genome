#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Test::Factory::DiskAllocation;

use Test::More tests => 2;

my $class = 'Genome::Disk::Command::Allocation::Transfer';
use_ok($class);

my $owner = UR::Value::Text->get('dummy-owner');
my $allocation = Genome::Test::Factory::DiskAllocation->setup_object(
    owner => $owner,
);

my $cmd = $class->create(
    allocations => [$allocation],
    target_group => $allocation->group,
    target_volume => $allocation->volume,
    remote_host => 'localhost', #this won't be a real test since we have no remote server
    remote_user => Genome::Sys->current_user->username,
    #remote_port => 22, #default
    authorization_key => 'dummy',
);
isa_ok($cmd, $class, 'instantiated class');

done_testing();
