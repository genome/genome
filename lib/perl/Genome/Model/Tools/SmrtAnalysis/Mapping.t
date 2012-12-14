#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

unless (Genome::Sys->username eq 'smrtanalysis') {
  plan skip_all => "this test is only runnable by user smrtanalysis"
}
plan tests => 3;

use above 'Genome';

use_ok('Genome::Model::Tools::SmrtAnalysis::Mapping');
my $ref_dir = '/gscmnt/pacbio/production/smrtanalysis/common/references/BAC_AC241402_3';
my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-Mapping';

my $input_fofn = $data_directory .'/input.fofn';
my $region_table = $data_directory .'/data/post_control_regions.fofn';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-Mapping-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 1,
);
my $control = Genome::Model::Tools::SmrtAnalysis::Mapping->create(
    job_directory => $tmp_dir,
    input_fofn => $input_fofn,
    post_control_fofn => $region_table,
    reference_directory => $ref_dir,
);
isa_ok($control,'Genome::Model::Tools::SmrtAnalysis::Mapping');
ok($control->execute,'Execute command '. $control->command_name);
