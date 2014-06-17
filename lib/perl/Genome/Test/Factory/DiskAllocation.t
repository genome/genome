#! /usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above "Genome";

use Test::More;

my $class = 'Genome::Test::Factory::DiskAllocation';
use_ok($class) or die;

ok(!eval{ $class->generate_obj; } , 'create test disk allocation');

class GenomeTest::Object { };
my $owner = GenomeTest::Object->create;
ok($owner, 'create owner');

my $disk_allocation = $class->generate_obj(owner => $owner);
ok($disk_allocation, 'create test disk allocation');
ok(-d $disk_allocation->absolute_path, 'disk_allocation absolute_path exists');

done_testing();
