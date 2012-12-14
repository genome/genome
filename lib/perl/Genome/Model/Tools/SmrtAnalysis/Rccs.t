#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

unless (Genome::Sys->username eq 'smrtanalysis') {
  plan skip_all => "this test is only runnable by user smrtanalysis"
}
plan tests => 3;

use_ok('Genome::Model::Tools::SmrtAnalysis::Rccs');

my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-Rccs';

my $input_fofn = $data_directory .'/input.fofn';
my $cmp_hdf5_file = $data_directory .'/data/aligned_reads.cmp.h5';
my $reference_directory = $ENV{SEYMOUR_HOME} .'/common/references/AML_121_amplicons';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-Rccs-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 0,
);
my $data_dir = $tmp_dir .'/data';
Genome::Sys->create_directory($data_dir);

my $rccs = Genome::Model::Tools::SmrtAnalysis::Rccs->create(
    input_fofn => $input_fofn,
    cmp_hdf5_file => $cmp_hdf5_file,
    job_directory => $tmp_dir,
    reference_directory => $reference_directory,
);
isa_ok($rccs,'Genome::Model::Tools::SmrtAnalysis::Rccs');
ok($rccs->execute,'Execute command '. $rccs->command_name);
