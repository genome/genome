#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

unless (Genome::Sys->username eq 'smrtanalysis') {
  plan skip_all => "this test is only runnable by user smrtanalysis"
}
plan tests => 3;

use_ok('Genome::Model::Tools::SmrtAnalysis::ConsensusReports');

my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-ConsensusReports';

my $reference_directory = '/gscmnt/pacbio/production/smrtanalysis/common/references/BAC_AC241402_3';
my $variants_gff = $data_directory .'/data/variants.gff';
my $alignment_summary_gff = $data_directory .'/data/alignment_summary.gff';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-ConsensusReports-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 1,
);
my $evi_cons = Genome::Model::Tools::SmrtAnalysis::ConsensusReports->create(
    reference_directory => $reference_directory,
    alignment_summary_gff_file => $alignment_summary_gff,
    variants_gff_file => $variants_gff,
    job_directory => $tmp_dir,
);
isa_ok($evi_cons,'Genome::Model::Tools::SmrtAnalysis::ConsensusReports');
ok($evi_cons->execute,'Execute command '. $evi_cons->command_name);

exit;
