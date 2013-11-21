#!/usr/bin/env genome-perl
use strict;
use warnings;

$ENV{UR_DBI_NO_COMMIT} = 1;

use Test::More;
use above "Genome";
use Carp::Always;

my $class = 'Genome::Config::Set';

use_ok($class);

my $config_set = Genome::Config::Set->create();
ok($config_set, 'it creates itself successfully');
isa_ok($config_set, $class, 'it returns the correct object type');
ok($config_set->allocation, 'it automatically creates an allocation if one is not supplied');
ok($config_set->path, 'it successfully delegates path to the underlying allocation');


my $allocation = Genome::Disk::Allocation->create(
    owner_id            => 'acoffman',
    disk_group_name     => $ENV{GENOME_DISK_GROUP_REFERENCES},
    allocation_path     => 'analysis_configuration/acoffman',
    owner_class_name    => 'Genome::Config::Set',
    kilobytes_requested => 25,
);
my $config_set_with_allocation = Genome::Config::Set->create(allocation => $allocation);

ok($config_set_with_allocation->allocation->id eq $allocation->id, 'it doesnt blow away an existing allocation if given one');

done_testing();

