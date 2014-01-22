#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Model::Tools::Somatic::TestHelpers qw( create_test_objects run_test );

my $TEST_DATA_VERSION = 9;

my $pkg = 'Genome::Model::Tools::Somatic::ProcessSomaticVariation';
use_ok($pkg);
my $main_dir = Genome::Utility::Test->data_dir_ok($pkg, $TEST_DATA_VERSION);

my $input_dir = File::Spec->join($main_dir, "input");

my $somatic_variation_build = create_test_objects($main_dir);

run_test(
    $pkg,
    $main_dir,
    somatic_variation_build => $somatic_variation_build,
    target_regions          => "$input_dir/target_regions.bed",
    igv_reference_name      => 'b37',
);

done_testing();
