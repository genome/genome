#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

if (Genome::Sys->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}
else {
    plan tests => 13;
}

use above 'Genome';

use_ok('Genome::Model::Tools::Vcf::Convert::Snv::VarscanSomatic');

my $test_dir = Genome::Config::get('test_inputs') . "/Genome-Model-Tools-Vcf-Convert-Snv-VarscanSomatic";
# V2 updated with various corrections
# V3 - remove VT INFO field in header and VT=SNP in the body
# V4 - add source in header
# V5 - Correct the AD and BQ fields number attribute
# V6 - Correct reference sequence specification and description of BQ and AD fields
# V7 - round VAQ to an integer
# V8 - VCF header format change
# V9 - Add DP4 and SS info, and VCF header change with center and tcgaversion
# v10- More TCGA-compliant vcf format
my $expected_base = "expected.v10";
my $expected_dir = "$test_dir/$expected_base";
my $expected_file = "$expected_dir/output.vcf";

my $output_file = Genome::Sys->create_temp_file_path;
my $input_file = "$test_dir/snvs.hq";

my %params = (
    input_file => $input_file, 
    output_file  => $output_file,
    aligned_reads_sample => "TEST",
    control_aligned_reads_sample => "OTHER_TEST",
    reference_sequence_build_id => 101947881,
);
my $command = Genome::Model::Tools::Vcf::Convert::Snv::VarscanSomatic->create(%params);

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

$command = Genome::Model::Tools::Vcf::Convert::Snv::VarscanSomatic->create(%params);
ok($command, 'Command created');
$rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file2, "TCGA output file created");

_run_if_ssl_new_enough( sub {
    compare_file($expected_file2, $output_file2);
} );

my $expected_file3 = $expected_dir . '/TCGA_output2.vcf';
my $output_file3 = Genome::Sys->create_temp_file_path;

$params{output_file} = $output_file3;
$params{aligned_reads_sample} = 'H_LS-BH-A0HQ-01A-11X-A049-09-1';
$params{control_aligned_reads_sample} = 'H_LS-AN-A046-10A-01X-A054-09-1';

$command = Genome::Model::Tools::Vcf::Convert::Snv::VarscanSomatic->create(%params);
ok($command, 'Command created');
$rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file3, "TCGA output file created");

_run_if_ssl_new_enough( sub {
    compare_file($expected_file3, $output_file3);
} );

done_testing();


sub _run_if_ssl_new_enough {
    my $sub = shift;

    use Net::SSLeay;
    SKIP: {
        if ($Net::SSLeay::VERSION < 1.74) {
            skip 'SSL is too old', 1;
        }
        $sub->();
    }
}

# The files will have a timestamp that will differ. Ignore this but check the rest.
sub compare_file {
    my ($expected_file, $out_file) = @_;
    my $expected = `cat $expected_file | grep -v fileDate`;
    my $output   = `zcat $out_file  | grep -v fileDate`;
    my $diff = Genome::Sys->diff_text_vs_text($output, $expected);
    ok(!$diff, 'output matched expected result') or diag("diff results:\n" . $diff);
    return 1;
}
