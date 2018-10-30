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

use_ok('Genome::Model::Tools::RefCov::WholeGenome');

my $tmp_dir = File::Temp::tempdir('BioSamtools-RefCov-'.Genome::Sys->username.'-XXXX',CLEANUP => 1, TMPDIR => 1);
my $data_dir = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-RefCov-WholeGenome';
my $expected_data_dir = $data_dir;

# INPUTS
my $alignment_file_path = $data_dir .'/NOTCH1_1k.bam';
my $regions_file = $data_dir .'/NOTCH1.bed';

# Excpected Outputs
my $expected_stats_file = $expected_data_dir .'/PDL_NOTCH1_1k_NOTCH1_STATS_expected.tsv';
my $expected_merged_stats_file = $expected_data_dir .'/PDL_NOTCH1_merged_stats_expected.tsv';
   
# OUTPUTS
my $merged_stats_file = $tmp_dir .'/NOTCH1_merged_stats.tsv';
my $stats_file = $tmp_dir .'/NOTCH1_1k_NOTCH1_STATS.tsv';

my $ref_build = Genome::Model::Build->get('101947881');
my $fasta = $ref_build->full_consensus_path('fa');

my $ref_cov = Genome::Model::Tools::RefCov::WholeGenome->create(
    stats_file => $stats_file,
    merged_stats_file => $merged_stats_file,
    merge_by => 'gene',
    alignment_file_path => $alignment_file_path,
    roi_file_path => $regions_file,
    #THIS TAKES TOO LONG FOR THIS TEST
    genome_normalized_coverage => 0,
    min_depth_filter => 1,
    reference_fasta => $fasta,
);
isa_ok($ref_cov,'Genome::Model::Tools::RefCov::WholeGenome');
ok($ref_cov->execute,'execute WholeGenome command '. $ref_cov->command_name);

ok(!compare($expected_stats_file,$ref_cov->stats_file),'expected stats file '. $expected_stats_file .' is identical to '. $ref_cov->stats_file);
ok(!compare($expected_merged_stats_file,$ref_cov->merged_stats_file),'expected merged stats file '. $expected_merged_stats_file .' is identical to '. $ref_cov->merged_stats_file);
