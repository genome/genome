#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use File::Compare;

use above 'Genome';

if ($] < 5.010) {
  plan skip_all => "this test is only runnable on perl 5.10+"
}
plan tests => 10;

use_ok('Genome::Model::Tools::RefCov::ExomeCapture');

my $tmp_dir = File::Temp::tempdir('BioSamtools-RefCov-'.Genome::Sys->username.'-XXXX',CLEANUP => 1, TMPDIR => 1);
my $data_dir = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-RefCov-ExomeCapture';
my $expected_data_dir = $data_dir;

my $alignment_file_path = $data_dir .'/test.bam';
my $regions_file = $data_dir .'/test.bed';

my $ref_build = Genome::Model::Build->get('101947881');
my $fasta = $ref_build->full_consensus_path('fa');

my $expected_stats_file = $expected_data_dir .'/PDL_test_STATS.tsv';
my $ref_cov = Genome::Model::Tools::RefCov::ExomeCapture->create(
    output_directory => $tmp_dir,
    alignment_file_path => $alignment_file_path,
    roi_file_path => $regions_file,
    evaluate_gc_content => 1,
    reference_fasta => $fasta,
);
isa_ok($ref_cov,'Genome::Model::Tools::RefCov::ExomeCapture');
ok($ref_cov->execute,'execute Standard command '. $ref_cov->command_name);

ok(!compare($expected_stats_file,$ref_cov->stats_file),'expected stats file '. $expected_stats_file .' is identical to '. $ref_cov->stats_file);
unlink($ref_cov->stats_file);

my $expected_q20_stats_file = $expected_data_dir .'/PDL_test_STATS-q20.tsv';
my $q20_ref_cov = Genome::Model::Tools::RefCov::ExomeCapture->create(
    output_directory => $tmp_dir,
    alignment_file_path => $alignment_file_path,
    roi_file_path => $regions_file,
    min_base_quality => 20,
    reference_fasta => $fasta,
);
isa_ok($q20_ref_cov,'Genome::Model::Tools::RefCov::ExomeCapture');
ok($q20_ref_cov->execute,'execute Standard command '. $q20_ref_cov->command_name);
ok(!compare($expected_q20_stats_file,$q20_ref_cov->stats_file),'expected stats file '. $expected_q20_stats_file .' is identical to '. $q20_ref_cov->stats_file);
unlink($q20_ref_cov->stats_file);

my $expected_q20_q1_stats_file = $expected_data_dir .'/PDL_test_STATS-q20-q1.tsv';
my $q20_q1_ref_cov = Genome::Model::Tools::RefCov::ExomeCapture->create(
    output_directory => $tmp_dir,
    alignment_file_path => $alignment_file_path,
    roi_file_path => $regions_file,
    min_base_quality => 20,
    min_mapping_quality => 1,
    reference_fasta => $fasta,
);
isa_ok($q20_q1_ref_cov,'Genome::Model::Tools::RefCov::ExomeCapture');
ok($q20_q1_ref_cov->execute,'execute Standard command '. $q20_q1_ref_cov->command_name);
ok(!compare($expected_q20_q1_stats_file,$q20_q1_ref_cov->stats_file),'expected stats file '. $expected_q20_q1_stats_file .' is identical to '. $q20_q1_ref_cov->stats_file);
