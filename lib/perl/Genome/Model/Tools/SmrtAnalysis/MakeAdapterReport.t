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

use_ok('Genome::Model::Tools::SmrtAnalysis::MakeAdapterReport');

my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-MakeAdapterReport';

my $hdf5_fofn = $data_directory .'/input.fofn';
my $expected_report_xml_file = $data_directory .'/results/filterReports_adapters.xml';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-MakeAdapterReport-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 1,
);
my $output_dir = $tmp_dir .'/results';
my $report_xml_file = $output_dir .'/filterReports_adapters.xml';
my $tool = Genome::Model::Tools::SmrtAnalysis::MakeAdapterReport->create(
    hdf5_fofn => $hdf5_fofn,
    report_xml_file => $report_xml_file,
    output_dir => $output_dir,
);
isa_ok($tool,'Genome::Model::Tools::SmrtAnalysis::MakeAdapterReport');
ok($tool->execute,'Execute command '. $tool->command_name);

is(compare($report_xml_file,$expected_report_xml_file),0,'Expected report xml file '. $expected_report_xml_file .' matches '. $report_xml_file);
