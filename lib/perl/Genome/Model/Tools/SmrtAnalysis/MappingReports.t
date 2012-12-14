#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

unless (Genome::Sys->username eq 'smrtanalysis') {
    plan skip_all => "this test is only runnable by user smrtanalysis"
}
plan tests => 3;

use above 'Genome';

use_ok('Genome::Model::Tools::SmrtAnalysis::MappingReports');

my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-MappingReports';

my $reference_directory = '/gscmnt/pacbio/production/smrtanalysis/common/references/BAC_AC241402_3';
my $alignment_summary_gff = $data_directory .'/data/alignment_summary.gff';
my $cmp_hdf5_file = $data_directory .'/data/aligned_reads.cmp.h5';
my $filtered_regions_fofn = $data_directory .'/data/filtered_regions.fofn';     

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-MappingReports-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 1,
);
my $filter = Genome::Model::Tools::SmrtAnalysis::MappingReports->create(
    reference_directory => $reference_directory,
    alignment_summary_gff => $alignment_summary_gff,
    job_directory => $tmp_dir,
    cmp_hdf5_file => $cmp_hdf5_file,
    filtered_regions_fofn => $filtered_regions_fofn,
);
isa_ok($filter,'Genome::Model::Tools::SmrtAnalysis::MappingReports');
ok($filter->execute,'Execute command '. $filter->command_name);
