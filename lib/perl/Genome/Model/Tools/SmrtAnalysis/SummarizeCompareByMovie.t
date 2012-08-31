#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;
use File::Compare;

unless (Genome::Sys->username eq 'smrtanalysis') {
  plan skip_all => "this test is only runnable by user smrtanalysis"
}
plan tests => 4;

use_ok('Genome::Model::Tools::SmrtAnalysis::SummarizeCompareByMovie');

my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-SummarizeCompareByMovie';

my $cmp_hdf5_file = $data_directory .'/data/control_reads.cmp.h5';
my $bas_hdf5_fofn = $data_directory .'/input.fofn';
my $expected_output_csv_file = $data_directory .'/data/control_results_by_movie.csv';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-SummarizeCompareByMovie-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 1,
);
my $output_dir = $tmp_dir .'/data';
Genome::Sys->create_directory($output_dir);
my $output_csv_file = $output_dir .'/control_results_by_movie.csv';
my $tool = Genome::Model::Tools::SmrtAnalysis::SummarizeCompareByMovie->create(
    cmp_hdf5_file => $cmp_hdf5_file,
    fofn => $bas_hdf5_fofn,
    output_csv_file => $output_csv_file,
    external => 1,
);
isa_ok($tool,'Genome::Model::Tools::SmrtAnalysis::SummarizeCompareByMovie');
ok($tool->execute,'Execute command '. $tool->command_name);
is(compare($output_csv_file,$expected_output_csv_file),0,'Expected output csv file '. $expected_output_csv_file .' matches the output csv file '. $output_csv_file);

exit;

