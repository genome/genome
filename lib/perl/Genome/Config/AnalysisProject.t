#!/usr/bin/env genome-perl
use strict;
use warnings;

$ENV{UR_DBI_NO_COMMIT} = 1;

use Test::More;
use above "Genome";
use Carp::Always;

my $class = 'Genome::Config::AnalysisProject';

use_ok($class);

my $test_obj = Genome::Config::AnalysisProject->create(
    created_by => Genome::Sys->username,
    name => 'Test Project',
    run_as => Genome::Sys->username,
);

isa_ok($test_obj, $class, 'It creates an object');

done_testing();


