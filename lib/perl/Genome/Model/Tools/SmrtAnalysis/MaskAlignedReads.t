#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

unless (Genome::Sys->username eq 'smrtanalysis') {
  plan skip_all => "this test is only runnable by user smrtanalysis"
}
plan tests => 3;

use_ok('Genome::Model::Tools::SmrtAnalysis::MaskAlignedReads');

my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-MaskAlignedReads';

my $cmp_hdf5_file = $data_directory .'/data/control_reads.cmp.h5';
my $region_table = $data_directory .'/data/filtered_regions.fofn';
my $expected_masked_table = $data_directory .'/data/post_control_regions.fofn';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-MaskAlignedReads-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 1,
);
my $output_dir = $tmp_dir .'/data';
Genome::Sys->create_directory($output_dir);
my $masked_table = $output_dir .'/post_control_regions.fofn';
my $tool = Genome::Model::Tools::SmrtAnalysis::MaskAlignedReads->create(
    cmp_hdf5_file => $cmp_hdf5_file,
    region_table => $region_table,
    masked_table => $masked_table,
);
isa_ok($tool,'Genome::Model::Tools::SmrtAnalysis::MaskAlignedReads');
ok($tool->execute,'Execute command '. $tool->command_name);
