#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use File::Compare;

use above 'Genome';

if ($] < 5.010) {
    plan skip_all => "this test is only runnable on perl 5.10+"
}
plan tests => 11;

use_ok('Genome::Model::Tools::BioSamtools');
use_ok('Genome::Model::Tools::BioSamtools::RefCov');

my $tmp_dir = File::Temp::tempdir('BioSamtools-RefCov-'.Genome::Sys->username.'-XXXX',DIR => "$ENV{GENOME_TEST_TEMP}",CLEANUP => 1);

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BioSamtools/RefCov';

my $bam_file = $data_dir .'/test.bam';
my $regions_file = $data_dir .'/test_regions_zero_based_start.bed';
my $expected_stats_file = $data_dir .'/test_test_regions_STATS-3.tsv';

my $ref_cov = Genome::Model::Tools::BioSamtools::RefCov->create(
    output_directory => $tmp_dir,
    bam_file => $bam_file,
    bed_file => $regions_file,
);
isa_ok($ref_cov,'Genome::Model::Tools::BioSamtools::RefCov');
ok($ref_cov->execute,'execute RefCov command '. $ref_cov->command_name);

ok(compare($expected_stats_file,$ref_cov->stats_file) == 0,'expected stats file '. $expected_stats_file .' is identical to '. $ref_cov->stats_file);
unlink($ref_cov->stats_file);

my $expected_q20_stats_file = $data_dir .'/test_test_regions_STATS-q20-2.tsv';
my $q20_ref_cov = Genome::Model::Tools::BioSamtools::RefCov->create(
    output_directory => $tmp_dir,
    bam_file => $bam_file,
    bed_file => $regions_file,
    min_base_quality => 20,
);
isa_ok($q20_ref_cov,'Genome::Model::Tools::BioSamtools::RefCov');
ok($q20_ref_cov->execute,'execute RefCov command '. $q20_ref_cov->command_name);
ok(compare($expected_q20_stats_file,$q20_ref_cov->stats_file) == 0,'expected stats file '. $expected_q20_stats_file .' is identical to '. $q20_ref_cov->stats_file);
unlink($q20_ref_cov->stats_file);

my $expected_q20_q1_stats_file = $data_dir .'/test_test_regions_STATS-q20-q1-2.tsv';
my $q20_q1_ref_cov = Genome::Model::Tools::BioSamtools::RefCov->create(
    output_directory => $tmp_dir,
    bam_file => $bam_file,
    bed_file => $regions_file,
    min_base_quality => 20,
    min_mapping_quality => 1,
);
isa_ok($q20_q1_ref_cov,'Genome::Model::Tools::BioSamtools::RefCov');
ok($q20_q1_ref_cov->execute,'execute RefCov command '. $q20_q1_ref_cov->command_name);
ok(compare($expected_q20_q1_stats_file,$q20_q1_ref_cov->stats_file) == 0,'expected stats file '. $expected_q20_q1_stats_file .' is identical to '. $q20_q1_ref_cov->stats_file);
