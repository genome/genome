#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;
use File::Compare;

unless (Genome::Sys->username eq 'smrtanalysis') {
  plan skip_all => "this test is only runnable by user smrtanalysis"
}
plan tests => 3;

use_ok('Genome::Model::Tools::SmrtAnalysis::CompareSequences');

my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-CompareSequences';
my $control_reference = $data_directory .'/Standard_v1';
my $bas_h5_fofn = $data_directory .'/input.fofn';
my $rgn_h5_fofn = $data_directory .'/output.fofn';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-CompareSequences-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 1,
);
my $data_dir = $tmp_dir .'/data';
unless (-d $data_dir) {
    Genome::Sys->create_directory($data_dir);
}
my $hdf5_file = $data_dir .'/control_reads.cmp.h5';
my $title_table = $data_dir .'/titleTable';
my $tool = Genome::Model::Tools::SmrtAnalysis::CompareSequences->create(
    algorithm => 'blasr',
    query => $bas_h5_fofn,
    region_table => $rgn_h5_fofn,
    target => $control_reference,
    split_subreads => 0,
    filter_adapter_only => 1,
    hdf5_file => $hdf5_file,
    hdf5_mode => 'w',
    algorithm_params => '-bestn 1 -titleTable '. $title_table,
    noise_data => '-77.27,0.08654,0.00121',
    min_z => 3,
    min_length => 50,
    min_accuracy => 0.75,
);
isa_ok($tool,'Genome::Model::Tools::SmrtAnalysis::CompareSequences');
ok($tool->execute,'Execute command '. $tool->command_name);
