#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

unless (Genome::Sys->username eq 'smrtanalysis') {
  plan skip_all => "this test is only runnable by user smrtanalysis"
}
plan tests => 3;

use_ok('Genome::Model::Tools::SmrtAnalysis::Resequencing');

my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-ConsensusStats';

my $alignment_summary_gff_file = $data_directory .'/data/alignment_summary.gff';
my $variants_gff_file = $data_directory .'/data/variants.gff';
my $cmp_hdf5_file = $data_directory .'/data/aligned_reads.cmp.h5';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-ConsensusStats-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 1,
);

my $reseq = Genome::Model::Tools::SmrtAnalysis::ConsensusStats->create(
    alignment_summary_gff_file => $alignment_summary_gff_file,
    variants_gff_file => $variants_gff_file,
    cmp_hdf5_file => $cmp_hdf5_file,
    output_gff_file => $tmp_dir .'/alignment_summary.gff',
);
isa_ok($reseq,'Genome::Model::Tools::SmrtAnalysis::ConsensusStats');
ok($reseq->execute,'Execute command '. $reseq->command_name);
