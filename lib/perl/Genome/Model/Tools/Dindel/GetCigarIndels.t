#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Model::Tools::Dindel::TestHelpers qw(
    get_test_dir
    get_ref_fasta
    compare_output_to_test_data
);

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

my $class = 'Genome::Model::Tools::Dindel::GetCigarIndels';
use_ok($class);

my $VERSION = 1; # Bump this each time test data changes

my $ref_fasta = get_ref_fasta();

my $test_dir = get_test_dir($class, $VERSION);
my $input_bam = File::Spec->join($test_dir, '134053361-region-limited.bam');
ok(-s $input_bam, "Found input Bam file");

my $output_directory = Genome::Sys->create_temp_directory();

my $cmd = $class->create(
    input_bam => $input_bam,
    ref_fasta => $ref_fasta,
    output_directory => $output_directory,
);
ok($cmd->execute(), "Successfully ran command");

compare_output_to_test_data($cmd->output_variants, $output_directory, $test_dir);
compare_output_to_test_data($cmd->output_libraries, $output_directory, $test_dir);

done_testing();
