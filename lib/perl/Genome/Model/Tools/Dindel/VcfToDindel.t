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

my $class = 'Genome::Model::Tools::Dindel::VcfToDindel';
use_ok($class);

my $VERSION = 0; # Bump this each time test data changes

my $ref_fasta = get_ref_fasta();

my $test_dir = get_test_dir($class, $VERSION);
my $input_vcf = File::Spec->join($test_dir, 'testdata-indel.vcf');
ok(-s $input_vcf, "Found input Vcf file") || die;

my $output_directory = Genome::Sys->create_temp_directory();

my $cmd = $class->create(
    input_vcf => $input_vcf,
    ref_fasta => $ref_fasta,
    output_directory => $output_directory,
);
ok($cmd->execute(), 'Successfully ran command') || die;

compare_output_to_test_data($cmd->output_dindel_file, $output_directory, $test_dir);

done_testing();
