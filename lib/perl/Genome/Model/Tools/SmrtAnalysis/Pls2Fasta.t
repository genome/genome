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

use_ok('Genome::Model::Tools::SmrtAnalysis::Pls2Fasta');

my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-Pls2Fasta';

my $bas_h5_fofn = $data_directory .'/input.fofn';
my $rgn_h5_fofn = $data_directory .'/output.fofn';
my $expected_fasta_file = $data_directory .'/filtered_subreads.fa';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-Pls2Fasta-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 1,
);
my $fasta_file = $tmp_dir .'/filtered_subreads.fa';
my $tool = Genome::Model::Tools::SmrtAnalysis::Pls2Fasta->create(
    hdf5_file => $bas_h5_fofn,
    region_table => $rgn_h5_fofn,
    fasta_file => $fasta_file,
    trim_by_region => 1,
);
isa_ok($tool,'Genome::Model::Tools::SmrtAnalysis::Pls2Fasta');
ok($tool->execute,'Execute command '. $tool->command_name);

is(compare($fasta_file,$expected_fasta_file),0,'Expected FASTA file '. $expected_fasta_file .' matches '. $fasta_file);

exit;

