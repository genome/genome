#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

my $class = 'Genome::Model::Tools::Dindel::VcfToDindel';
use_ok($class);

my $VERSION = 0; # Bump this each time test data changes

my $test_dir = File::Spec->join(Genome::Utility::Test->data_dir($class), "v$VERSION");
diag "Test data located at $test_dir\n";

my $input_vcf = File::Spec->join($test_dir, 'testdata-indel.vcf');
ok(-s $input_vcf, "Found input Vcf file") || die;

my $reference_sequence_build = Genome::Model::Build->get(106942997);
my $ref_fasta = $reference_sequence_build->full_consensus_path('fa');
diag "Reference sequence found at: $ref_fasta\n";
ok(-s $ref_fasta, "Found Reference Sequence fasta") || die;

my $output_directory = Genome::Sys->create_temp_directory();

my $cmd = $class->create(
    input_vcf => $input_vcf,
    ref_fasta => $ref_fasta,
    output_directory => $output_directory,
);
ok($cmd->execute(), 'Successfully ran command');

my $expected_output_file = $cmd->output_dindel_file;
$expected_output_file =~ s/$output_directory/$test_dir/;
compare_ok($expected_output_file, $cmd->output_dindel_file, 'Output matches expected output');

done_testing();
