#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use File::Compare;

use above 'Genome';

if ($] < 5.010) {
  plan skip_all => "this test is only runnable on perl 5.10+"
}
plan tests => 5;

use_ok('Genome::Model::Tools::BioSamtools');
use_ok('Genome::Model::Tools::BioSamtools::CompareAlignmentSummaries');


my $tmp_dir = File::Temp::tempdir('BioSamtools-CompareAlignmentSummaries-'.Genome::Sys->username.'-XXXX',CLEANUP => 1, TMPDIR => 1);

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BioSamtools/AlignmentSummary';
my $input_file_1 = $data_dir .'/alignment_summary.tsv';
my $input_file_2 = $data_dir .'/alignment_summary_2.tsv';
my $output_file = $tmp_dir .'/merged_alignment_summary.tsv';
my $expected_output_file = $data_dir .'/merged_alignment_summary-4.tsv';

my $mas = Genome::Model::Tools::BioSamtools::CompareAlignmentSummaries->create(
    output_file => $output_file,
    input_files => [$input_file_1,$input_file_2],
);
isa_ok($mas,'Genome::Model::Tools::BioSamtools::CompareAlignmentSummaries');
ok($mas->execute,'execute AlignmentSummary command '. $mas->command_name);

ok(!compare($expected_output_file,$mas->output_file),'expected output file '. $expected_output_file .' is identical to '. $mas->output_file);
