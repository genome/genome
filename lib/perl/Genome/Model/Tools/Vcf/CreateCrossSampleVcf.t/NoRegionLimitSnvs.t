#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Model::Tools::Vcf::CreateCrossSampleVcf::TestHelpers qw(
    test_cmd
);

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
};

Genome::Config::set_env('workflow_builder_backend', 'inline');

my $VERSION = 'no-region-limit-snvs-3';
my $use_mg = 0;
my $no_region_limiting = 1;
test_cmd('snvs', $VERSION, $use_mg, $no_region_limiting);

done_testing();
