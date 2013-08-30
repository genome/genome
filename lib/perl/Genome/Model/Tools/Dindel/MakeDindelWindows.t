#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Model::Tools::Dindel::TestHelpers qw(
    get_test_dir
    compare_output_to_test_data
);

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

my $class = 'Genome::Model::Tools::Dindel::MakeDindelWindows';
use_ok($class);

my $VERSION = 0; # Bump this each time tests data changes

my $test_dir = get_test_dir($class, $VERSION);

my $input_dindel_file = File::Spec->join($test_dir, 'input.dindel');
ok(-s $input_dindel_file, 'Found input dindel file') || die;

my $output_directory = Genome::Sys->create_temp_directory();

my $cmd = $class->create(
    input_dindel_file => $input_dindel_file,
    output_directory => $output_directory,
);
ok($cmd->execute(), "Successfully ran command") || die;

compare_output_to_test_data($cmd->window_file, $output_directory, $test_dir);

done_testing();
