#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;
use File::Compare;

unless (Genome::Sys->username eq 'smrtanalysis') {
  plan skip_all => "this test is only runnable by user smrtanalysis"
}
plan tests => 13;

use_ok('Genome::Model::Tools::SmrtAnalysis::FilterPlsH5');

my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-FilterPlsH5';
my @bas_h5_files = glob($data_directory .'/*.bas.h5');

my $expected_fofn_file_path = $data_directory .'/input.fofn';
my $expected_output_summary = $data_directory .'/summary.csv';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-FilterPlsH5-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 0,
);
my $fofn_file_path = $tmp_dir .'/input.fofn';

my $fofn_fh = Genome::Sys->open_file_for_writing($fofn_file_path);
for my $bas_h5_file (@bas_h5_files) {
    print $fofn_fh $bas_h5_file ."\n";
}
$fofn_fh->close;

is(compare($fofn_file_path,$expected_fofn_file_path),0,'Expected input.fofn file '. $expected_fofn_file_path .' matches '. $fofn_file_path);

my $output_dir = $tmp_dir .'/filtered_regions';
my $output_summary = $tmp_dir .'/summary.csv';
my $output_fofn = $tmp_dir .'/output.fofn';
my $tool = Genome::Model::Tools::SmrtAnalysis::FilterPlsH5->create(
    input_fofn => $fofn_file_path,
    min_read_length => 50,
    min_read_score => 0.75,
    output_directory => $output_dir,
    output_summary => $output_summary,
    output_fofn => $output_fofn,
);
isa_ok($tool,'Genome::Model::Tools::SmrtAnalysis::FilterPlsH5');
ok($tool->execute,'Execute command '. $tool->command_name);

is(compare($output_summary,$expected_output_summary),0,'Expected summary.csv file '. $expected_output_summary .' matches '. $output_summary);

ok(-s $output_fofn,'Output fofn'. $output_fofn .' has size.');
my $fh = Genome::Sys->open_file_for_reading($output_fofn);

my @rgn_h5_files = glob($output_dir .'/*rgn.h5');
ok((scalar(@bas_h5_files) == scalar(@rgn_h5_files)),'Found 1 region HDF5 for each base HDF5');
for my $rgn_h5_file (@rgn_h5_files) {
    my $line = $fh->getline;
    chomp($line);
    ok(($line eq $rgn_h5_file),'Output fofn line '. $line .' matches expected '. $rgn_h5_file);
    ok(-s $rgn_h5_file,'Output region HDF5 file has size.');
}
$fh->close;

# The second test is designed to use the read_white_list feature and RCCS reads which may require more stringent params
my $rccs_data_directory = $data_directory .'/RCCS';
my $rccs_input_fofn = $rccs_data_directory .'/input.fofn';
my $rccs_output_directory = $tmp_dir .'/RCCS/data';
Genome::Sys->create_directory($rccs_output_directory);
my $rccs_tool = Genome::Model::Tools::SmrtAnalysis::FilterPlsH5->create(
    input_fofn => $rccs_input_fofn,
    min_read_length => 50,
    min_read_score => 0.75,
    base_output_directory => $rccs_output_directory,
);
isa_ok($rccs_tool,'Genome::Model::Tools::SmrtAnalysis::FilterPlsH5');
ok($rccs_tool->execute,'Execute command '. $rccs_tool->command_name);

my $rccs_read_white_list_output_directory = $tmp_dir .'/white_list_data';
Genome::Sys->create_directory($rccs_read_white_list_output_directory);
my $rccs_read_white_list = $rccs_data_directory .'/read_white_list.txt';
my $rccs_read_white_list_tool = Genome::Model::Tools::SmrtAnalysis::FilterPlsH5->create(
    input_fofn => $rccs_input_fofn,
    min_read_length => 50,
    min_read_score => 0.75,
    read_white_list => $rccs_read_white_list,
    base_output_directory => $rccs_read_white_list_output_directory,
);
isa_ok($rccs_read_white_list_tool,'Genome::Model::Tools::SmrtAnalysis::FilterPlsH5');
ok($rccs_read_white_list_tool->execute,'Execute command '. $rccs_read_white_list_tool->command_name);

exit;

