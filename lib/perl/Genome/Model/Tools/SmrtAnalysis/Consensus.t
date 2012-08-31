#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

unless (Genome::Sys->username eq 'smrtanalysis') {
  plan skip_all => "this test is only runnable by user smrtanalysis"
}
plan tests => 3;

use_ok('Genome::Model::Tools::SmrtAnalysis::Consensus');

my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-Consensus';

my $original_cmp_hdf5_file = $data_directory .'/data/aligned_reads.cmp.h5';
my $original_alignment_summary_gff = $data_directory .'/data/alignment_summary.gff';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-Consensus-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 0,
);

my $data_dir = $tmp_dir .'/data';
Genome::Sys->create_directory($data_dir);
my $alignment_summary_gff = $data_dir .'/alignment_summary.gff';
File::Copy::copy($original_alignment_summary_gff,$alignment_summary_gff) || die;

my $cmp_hdf5_file = $data_dir .'/aligned_reads.cmp.h5';
File::Copy::copy($original_cmp_hdf5_file,$cmp_hdf5_file) || die;

my $evi_cons = Genome::Model::Tools::SmrtAnalysis::Consensus->create(
    cmp_hdf5_file => $cmp_hdf5_file,
    alignment_summary_gff => $alignment_summary_gff,
    job_directory => $tmp_dir,
    nproc => 1,
);
isa_ok($evi_cons,'Genome::Model::Tools::SmrtAnalysis::Consensus');
ok($evi_cons->execute,'Execute command '. $evi_cons->command_name);

exit;
