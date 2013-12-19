#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 5;

use_ok('Genome::Model::Tools::Vcf::Backfill');

my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Vcf-Backfill";
my $input_dir = "$test_dir/input.v3";
my $expected_base = "expected.v5";
my $expected_dir = "$test_dir/$expected_base";
my $expected_file = "$expected_dir/output.vcf";

my $output_file = Genome::Sys->create_temp_file_path;
my $input_vcf= "$input_dir/snvs.H_LH-ED0274-25330081.region_limited.vcf.gz";
my $input_pileup = "$input_dir/H_LH-ED0274-25330081.for_snvs.pileup";
my $input_merged = "$input_dir/merged_positions.bed.gz";

my $command = Genome::Model::Tools::Vcf::Backfill->create( vcf_file => $input_vcf, 
                                                           pileup_file  => $input_pileup,
                                                           output_file  => $output_file,
                                                           merged_positions_file => $input_merged,
                                                       );

ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file created");

# The files will have a timestamp that will differ. Ignore this but check the rest.
my $diff = Genome::Sys->diff_file_vs_file($output_file, $expected_file);
#print "$output_file $expected_file\n"; sleep 10000;
ok(!$diff, 'output matched expected result')
    or diag("diff results:\n" . $diff);

