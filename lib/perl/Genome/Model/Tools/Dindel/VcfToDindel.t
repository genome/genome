#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 5;
use Genome::Utility::Test;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

my $class = 'Genome::Model::Tools::Dindel::VcfToDindel';
use_ok($class);

my $VERSION = 0; # Bump this each time test data changes

my $test_dir = File::Spec->join(Genome::Utility::Test->data_dir($class), "v$VERSION");
diag "Test data located at $test_dir\n";

my $input_vcf = File::Spec->join($test_dir, 'testdata-indel.vcf');
ok(-s $input_vcf, "Found input Vcf file");

my $reference_sequence_build = Genome::Model::Build->get(106942997);
my $ref_fasta = $reference_sequence_build->full_consensus_path('fa');
diag "Reference sequence found at: $ref_fasta\n";
ok(-s $ref_fasta, "Found Reference Sequence fasta");

my $output_file = Genome::Sys->create_temp_file_path();

my $cmd = $class->create(
    input_vcf => $input_vcf,
    output_dindel_file => $output_file,
    ref_fasta => $ref_fasta,
);
ok($cmd->execute(), 'Successfully ran command');

my $expected_output_file = File::Spec->join($test_dir, 'output_dindel_file.out');
my $diff = `diff $output_file $expected_output_file`;
ok(!$diff, 'Output files are identical');
