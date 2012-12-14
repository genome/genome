#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

use above 'Genome';

unless (Genome::Sys->username eq 'smrtanalysis') {
  plan skip_all => "this test is only runnable by user smrtanalysis"
}
plan tests => 3;

use_ok('Genome::Model::Tools::SmrtAnalysis::FilterReports');

my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-FilterReports';
my $input_fofn = $data_directory .'/input.fofn';
my $filtered_summary_file = $data_directory .'/data/filtered_summary.csv';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-FilterReports-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 1,
);
my $filter = Genome::Model::Tools::SmrtAnalysis::FilterReports->create(
    job_directory => $tmp_dir,
    input_fofn => $input_fofn,
    filtered_summary_file => $filtered_summary_file,
);
isa_ok($filter,'Genome::Model::Tools::SmrtAnalysis::FilterReports');
ok($filter->execute,'Execute command '. $filter->command_name);
