#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Utility::Vcf ('diff_vcf_file_vs_file');
use Test::More tests => 5;

use_ok('Genome::Model::Tools::Beagle::PhaseVcf');

my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Beagle-PhaseVcf";
my $input_dir = "$test_dir/input.v1";
my $input_vcf = "$input_dir/input.vcf.gz";
my $input_ped = "$input_dir/input.ped";
my $expected_dir = "$test_dir/expected.v1";
my $expected_vcf = "$expected_dir/expected.vcf.gz";
my $output_file_base = Genome::Sys->create_temp_file_path;
my $output_vcf = $output_file_base.".vcf.gz";
my $chrom = "1";

my $command = Genome::Model::Tools::Beagle::PhaseVcf->create( vcf_file => $input_vcf, 
    output_file  => $output_file_base,
    ped_file => $input_ped,
    chromosome => $chrom,
);

ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_vcf, "output file created");

my $diff = diff_vcf_file_vs_file($output_vcf, $expected_vcf);
ok(!$diff, "output $output_vcf matched expected result $expected_vcf")
    or diag("diff results:\n" . $diff);


