#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use File::Compare;
use File::Temp;

use above 'Genome';

if ($] < 5.010) {
  plan skip_all => "this test is only runnable on perl 5.10+"
}
plan tests => 4;

use_ok('Genome::Model::Tools::BioSamtools::MultiAlignmentStats');

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BioSamtools-MultiAlignmentStats';
my $expected_data_dir = $data_dir;

my $alignment_file_path = $data_dir .'/test.bam';
my $expected_bam_file = $data_dir .'/unique_test.bam';
my $expected_stats_file = $expected_data_dir .'/stats.tsv';

my $tmp_tsv = $tmp_dir .'/stats.tsv';
my $tmp_unique_bam = $tmp_dir .'/unique_test.bam';
my $stats = Genome::Model::Tools::BioSamtools::MultiAlignmentStats->create(
    output_stats_tsv => $tmp_tsv,
    aligned_bam_file => $alignment_file_path,
    unique_bam_file => $tmp_unique_bam,
);
isa_ok($stats,'Genome::Model::Tools::BioSamtools::MultiAlignmentStats');
ok($stats->execute,'execute Standard command '. $stats->command_name);
ok(!compare($expected_stats_file,$stats->output_stats_tsv),'expected stats file '. $expected_stats_file .' is identical to '. $stats->output_stats_tsv);
