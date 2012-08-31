#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use File::Compare;

use above 'Genome';

if ($] < 5.010) {
  plan skip_all => "this test is only runnable on perl 5.10+"
}
plan tests => 6;

use_ok('Genome::Model::Tools::BioSamtools');
use_ok('Genome::Model::Tools::BioSamtools::CleanBam');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BioSamtools-CleanBam';

my $input_bam_file = $data_dir .'/dirty.bam';
my $output_bam_file = Genome::Sys->create_temp_file_path('clean.bam');  
my $expected_bam_file = $data_dir .'/clean.bam';
my $summary_output_file = Genome::Sys->create_temp_file_path('summary.txt');    
my $expected_summary_output_file = $data_dir .'/summary.txt';
my $cb = Genome::Model::Tools::BioSamtools::CleanBam->create(
   input_bam_file => $input_bam_file,   
   output_bam_file => $output_bam_file,
   summary_output_file => $summary_output_file,
);
isa_ok($cb,'Genome::Model::Tools::BioSamtools::CleanBam');
ok($cb->execute,'execute CleanBam command '. $cb->command_name);
ok(!compare($expected_bam_file,$output_bam_file),'expected output BAM file '. $expected_bam_file .' is identical to '. $output_bam_file);
ok(!compare($expected_summary_output_file,$summary_output_file),'expected output BAM file '. $expected_summary_output_file .' is identical to '. $summary_output_file);
exit;
