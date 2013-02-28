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

# run_workflow_lsf uses this environment variable to allow it to run inline.
$ENV{NO_LSF} = 1;

use_ok('Genome::Model::Tools::BioSamtools');
use_ok('Genome::Model::Tools::BioSamtools::RefCov');
use_ok('Genome::Model::Tools::BioSamtools::ParallelRefCov');

my $tmp_dir = File::Temp::tempdir('BioSamtools-RefCov-'.Genome::Sys->username.'-XXXX',CLEANUP => 1, TMPDIR => 1);

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BioSamtools/RefCov';

my $bam_file = $data_dir .'/test.bam';
my $bed_file = $data_dir .'/test_regions_zero_based_start.bed';
my $expected_stats_file = $data_dir .'/test_test_regions_STATS-3.tsv';

my $ref_cov = Genome::Model::Tools::BioSamtools::ParallelRefCov->create(
    output_directory => $tmp_dir,
    bam_file => $bam_file,
    bed_file => $bed_file,
    regions => 5,
);
isa_ok($ref_cov,'Genome::Model::Tools::BioSamtools::ParallelRefCov');
ok($ref_cov->execute,'execute RefCov command '. $ref_cov->command_name);

ok(!compare($expected_stats_file,$ref_cov->stats_file),'expected stats file '. $expected_stats_file .' is identical to '. $ref_cov->stats_file);
