#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}
else {
    plan tests => 9;
}

use above 'Genome';

use_ok('Genome::Model::Tools::Vcf::Convert::Snv::Sniper');

my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Vcf-Convert-Snv-Sniper";
# V2 updated with various corrections
# V3 - remove VT INFO field in header and VT=SNP in the body
# V4 - add source in header
# V5 - Correct the AD and BQ fields number attribute
# V6 - Correct reference sequence specification and description of BQ and AD fields
# V7 - VCF header change
# v8 - VCF header change: center and tcgaversion
# v9 - VCF header change: TCGA-compliant
my $expected_base = "expected.v9";
my $expected_dir = "$test_dir/$expected_base";
my $expected_file = "$expected_dir/output.vcf";

my $output_file = Genome::Sys->create_temp_file_path;
my $input_file = "$test_dir/snvs.hq";

my %params = (
    input_file  => $input_file, 
    output_file => $output_file,
    aligned_reads_sample => "TUMOR_SAMPLE_123",
    control_aligned_reads_sample => "CONTROL_SAMPLE_123",
    reference_sequence_build_id  => 101947881,
);
my $command = Genome::Model::Tools::Vcf::Convert::Snv::Sniper->create(%params);

ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file created");

# The files will have a timestamp that will differ. Ignore this but check the rest.
my $expected = `cat $expected_file | grep -v fileDate`;
my $output = `zcat $output_file | grep -v fileDate`;
my $diff = Genome::Sys->diff_text_vs_text($output, $expected);
ok(!$diff, 'output matched expected result')
    or diag("diff results:\n" . $diff);


#This is to test new version of Sniper that can produce vcf as raw output
$expected_file = "$expected_dir/raw_vcf_output.vcf";
$output_file   = Genome::Sys->create_temp_file_path;

$params{input_file}  = "$test_dir/vcf-test.v1/snvs.hq";
$params{output_file} = $output_file;
$params{reference_sequence_build_id} = 106942997;

$command = Genome::Model::Tools::Vcf::Convert::Snv::Sniper->create(%params);

ok($command, 'Command created');
$rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file $output_file created");

# The files will have a timestamp that will differ. Ignore this but check the rest.
$expected = `cat $expected_file | grep -v fileDate`;
$output   = `zcat $output_file | grep -v fileDate`;
$diff = Genome::Sys->diff_text_vs_text($output, $expected);
ok(!$diff, 'output matched expected result')
    or diag("diff results:\n" . $diff);

