#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

unless (Genome::Sys->username eq 'smrtanalysis') {
  plan skip_all => "this test is only runnable by user smrtanalysis"
}
plan tests => 3;

use_ok('Genome::Model::Tools::SmrtAnalysis::EviConsWrapper');

my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-EviConsWrapper';

my $cmp_hdf5_file = $data_directory .'/data/aligned_reads.cmp.h5';
my $evi_cons_params = '--fastMode --baseMap=3,4,2,1 --runDecode --decodeFile=/gscmnt/pacbio/production/smrtanalysis/analysis/etc/defaultDecode.params --nproc=1 --hdf5Reference=ref000001 --subAlignment --refStart=0 --refEnd=99999';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-EviConsWrapper-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 1,
);
my $evi_cons = Genome::Model::Tools::SmrtAnalysis::EviConsWrapper->create(
    cmp_hdf5_file => $cmp_hdf5_file,
    evi_cons_params => $evi_cons_params,
    base_output_directory => $tmp_dir, 
);
isa_ok($evi_cons,'Genome::Model::Tools::SmrtAnalysis::EviConsWrapper');
ok($evi_cons->execute,'Execute command '. $evi_cons->command_name);

exit;
