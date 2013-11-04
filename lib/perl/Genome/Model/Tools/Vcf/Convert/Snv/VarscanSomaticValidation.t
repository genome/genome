#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}
else {
    plan tests => 13;
}

use above 'Genome';

use_ok('Genome::Model::Tools::Vcf::Convert::Snv::VarscanSomaticValidation');

my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Vcf-Convert-Snv-VarscanSomaticValidation";
#this is just a copy of the varscan somatic tests, with the header changed to include validation
#since everything important is inherited from that module

my $expected_base = "expected.v1";
my $expected_dir = "$test_dir/$expected_base";
my $expected_file = "$expected_dir/output.vcf";

my $output_file = Genome::Sys->create_temp_file_path;
my $input_file = "$test_dir/snvs.hq";
print STDERR $input_file . "\n";
my %params = (
    input_file => $input_file, 
    output_file  => $output_file,
    aligned_reads_sample => "TEST",
    control_aligned_reads_sample => "OTHER_TEST",
    reference_sequence_build_id => 101947881,
);
my $command = Genome::Model::Tools::Vcf::Convert::Snv::VarscanSomaticValidation->create(%params);

ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file created");
compare_file($expected_file, $output_file);

#test TCGA vcf
my $expected_file2 = $expected_dir . '/TCGA_output.vcf';
my $output_file2 = Genome::Sys->create_temp_file_path;

$params{output_file} = $output_file2;
$params{aligned_reads_sample} = 'TCGA-DS-A1OC-01A-11D-A14W-08';
$params{control_aligned_reads_sample} = 'TCGA-DS-A1OC-10A-01D-A14W-08';

$command = Genome::Model::Tools::Vcf::Convert::Snv::VarscanSomaticValidation->create(%params);
ok($command, 'Command created');
$rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file2, "TCGA output file created");

compare_file($expected_file2, $output_file2);

my $expected_file3 = $expected_dir . '/TCGA_output2.vcf';
my $output_file3 = Genome::Sys->create_temp_file_path;

$params{output_file} = $output_file3;
$params{aligned_reads_sample} = 'H_LS-BH-A0HQ-01A-11X-A049-09-1';
$params{control_aligned_reads_sample} = 'H_LS-AN-A046-10A-01X-A054-09-1';

$command = Genome::Model::Tools::Vcf::Convert::Snv::VarscanSomaticValidation->create(%params);
ok($command, 'Command created');
$rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file3, "TCGA output file created");
compare_file($expected_file3, $output_file3);

done_testing();


# The files will have a timestamp that will differ. Ignore this but check the rest.
sub compare_file {
    my ($expected_file, $out_file) = @_;
    my $expected = `cat $expected_file | grep -v fileDate`;
    my $output   = `zcat $out_file  | grep -v fileDate`;
    my $diff = Genome::Sys->diff_text_vs_text($output, $expected);
    ok(!$diff, 'output matched expected result') or diag("diff results:\n" . $diff);
    return 1;
}
