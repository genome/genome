#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use File::Temp qw/ tempdir /;
use File::Compare;

if ($] < 5.010) {
    plan skip_all => "this test is only runnable on perl 5.10+"
}
plan tests => 5;

my $test_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BioSamtools/Breakdown';
my $bam_file = $test_data_dir .'/breakdown.bam';

my $expected_tsv = $test_data_dir .'/breakdown-64.tsv';
my $base_output_dir = tempdir('BREAKDOWN_XXXXX',DIR => '/tmp',CLEANUP=> 1);

my $cmd = Genome::Model::Tools::BioSamtools::Breakdown->create(
    output_file => $base_output_dir .'/breakdown.tsv',
    bam_file => $bam_file,
);

isa_ok($cmd,'Genome::Model::Tools::BioSamtools::Breakdown');
ok($cmd->execute,'execute breakdown command '. $cmd->command_name);
ok(-f $cmd->output_file,'found output tsv file '. $cmd->output_file);
ok(-s $cmd->output_file,'output tsv file '. $cmd->output_file .' has size');
ok(!compare($cmd->output_file,$expected_tsv),'output tsv '. $cmd->output_file .' matches expected '. $expected_tsv);
