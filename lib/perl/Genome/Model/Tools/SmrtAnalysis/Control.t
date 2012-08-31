#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

use above 'Genome';

unless (Genome::Sys->username eq 'smrtanalysis') {
  plan skip_all => "this test is only runnable by user smrtanalysis"
}
plan tests => 3;

use_ok('Genome::Model::Tools::SmrtAnalysis::Control');
my $standard_dir = '/gscmnt/pacbio/production/smrtanalysis/common/references/Standard_v1';
my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-Control';

my $input_fofn = $data_directory .'/input.fofn';
my $region_table = $data_directory .'/data/filtered_regions.fofn';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-Control-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 1,
);
my $control = Genome::Model::Tools::SmrtAnalysis::Control->create(
    job_directory => $tmp_dir,
    input_fofn => $input_fofn,
    filtered_fofn => $region_table,
    control_reference_directory => $standard_dir,
);
isa_ok($control,'Genome::Model::Tools::SmrtAnalysis::Control');
ok($control->execute,'Execute command '. $control->command_name);

exit;
