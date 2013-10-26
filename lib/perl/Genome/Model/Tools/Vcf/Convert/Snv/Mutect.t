#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

use above 'Genome';

use_ok('Genome::Model::Tools::Vcf::Convert::Snv::Mutect');

#my $test_dir = "/gscmnt/sata831/info/medseq/dlarson/mutect_testing/Genome-Model-Tools-Vcf-Convert-Snv-Mutect";
my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Vcf-Convert-Snv-Mutect";

my $expected_base = "expected.v1";
my $expected_dir = "$test_dir/$expected_base";
my $expected_file = "$expected_dir/expected.vcf.gz";

my $input_file = "$test_dir/snvs.hq";
my $output_file = Genome::Sys->create_temp_file_path;
my %params = ( aligned_reads_sample => 'TestTumor',
               control_aligned_reads_sample => 'TestNormal',
               input_file => $input_file,
               output_file => $output_file,
               reference_sequence_build_id => 106942997,
           );


my $command = Genome::Model::Tools::Vcf::Convert::Snv::Mutect->create(%params);

ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file created");

# The files will have a timestamp that will differ. Ignore this but check the rest.
my $expected = `zcat $expected_file | grep -v fileDate`;
my $output = `zcat $output_file | grep -v fileDate`;
my $diff = Genome::Sys->diff_text_vs_text($output, $expected);
ok(!$diff, 'output matched expected result')
    or diag("diff results:\n" . $diff);
done_testing()
#perl -I ~/src/genome/lib/perl/ -S gmt vcf convert snv mutect --aligned-reads-sample 'TestTumor' --control-aligned-reads-sample 'TestNormal' --input-file snvs.hq --output-file expected.vcf.gz --reference-sequence-build-id 106942997
