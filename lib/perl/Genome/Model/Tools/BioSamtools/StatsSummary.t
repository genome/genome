#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use File::Compare;

use above 'Genome';

if ($] < 5.010) {
    plan skip_all => "this test is only runnable on perl 5.10+"
}
plan tests => 5;

use_ok('Genome::Model::Tools::BioSamtools');
use_ok('Genome::Model::Tools::BioSamtools::StatsSummary');

my $tmp_dir = File::Temp::tempdir('BioSamtools-StatsSummary-'.Genome::Sys->username.'-XXXX',CLEANUP => 1, TMPDIR => 1);
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BioSamtools/StatsSummary';

my $stats_file = $data_dir .'/test.stats';
my $expected_file = $data_dir .'/test_stats_summary-4.tsv';

my $as = Genome::Model::Tools::BioSamtools::StatsSummary->create(
    stats_file => $stats_file,
    output_file => $tmp_dir .'/test_stats_summary.tsv',
);
isa_ok($as,'Genome::Model::Tools::BioSamtools::StatsSummary');
ok($as->execute,'execute StatsSummary command '. $as->command_name);
ok(!(compare($expected_file,$as->output_file)),'expected file '. $expected_file .' matched the output file '. $as->output_file);
