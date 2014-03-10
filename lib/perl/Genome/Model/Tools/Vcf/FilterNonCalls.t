#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Utility::Vcf "diff_vcf_file_vs_file";
use Test::More;

my $class = 'Genome::Model::Tools::Vcf::FilterNonCalls';
use_ok($class);

my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Vcf-FilterNonCalls";
my $expected_base = "expected";
my $expected_dir = "$test_dir/$expected_base";
my $expected_file = "$expected_dir/expected.vcf.gz";

my $output_file = Genome::Sys->create_temp_file_path;
my $input_file = "$test_dir/input/snvs.vcf.gz";
ok(-s $input_file, "input file exists $input_file");

my $cmd = $class->create(input_file => $input_file, output_file => $output_file);

ok($cmd, 'Command created');
my $rv = $cmd->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file created $output_file");
ok(-s $expected_file, "expected file exists $expected_file");

my $diff = diff_vcf_file_vs_file($output_file, $expected_file);
ok(!$diff, 'Output matches expected result') 
    or diag("diff results:\n" . $diff);

done_testing();
