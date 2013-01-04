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

my $data_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SmrtAnalysis-Resequencing';

my $reference_directory = $data_directory .'/references/lambda';
my $control_reference_directory = $data_directory .'/references/Standard_v1';
my $input_fofn = $data_directory .'/input.fofn';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-SmrtAnalysis-Resequencing-XXXXXX',
    DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
    CLEANUP => 1,
);

my $reseq = Genome::Model::Tools::SmrtAnalysis::Resequencing->create(
    reference_directory => $reference_directory,
    control_reference_directory => $control_reference_directory,
    input_fofn =>  $input_fofn,
    job_directory => $tmp_dir,
);
isa_ok($reseq,'Genome::Model::Tools::SmrtAnalysis::Resequencing');
ok($reseq->execute,'Execute command '. $reseq->command_name);
