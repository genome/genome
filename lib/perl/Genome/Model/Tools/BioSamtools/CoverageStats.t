#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

use above 'Genome';

if ($] < 5.010) {
    plan skip_all => "this test is only runnable on perl 5.10+"
}
plan tests => 6;

Genome::Config::set_env('workflow_builder_backend', 'inline');

map { print STDERR $_ . " => " . $ENV{$_} . "\n" } keys %ENV;

use_ok('Genome::Model::Tools::BioSamtools::CoverageStats');

my $tmp_dir = File::Temp::tempdir('BioSamtools-CoverageStats-'.Genome::Sys->username.'-XXXX',CLEANUP => 1, TMPDIR => 1);

my $data_dir = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-BioSamtools/RefCov';

my $bam_file = $data_dir .'/test.bam';
my $regions_file = $data_dir .'/test_regions.bed';
my @wingspans = qw/0 200/;
my @minimum_depths = qw/1 10 20/;
my $stats = Genome::Model::Tools::BioSamtools::CoverageStats->create(
    output_directory => $tmp_dir,
    bam_file => $bam_file,
    bed_file => $regions_file,
    wingspan_values => join(',',@wingspans),
    minimum_depths => join(',',@minimum_depths),
    minimum_base_quality => 20,
    minimum_mapping_quality => 1,
);
isa_ok($stats,'Genome::Model::Tools::BioSamtools::CoverageStats');
ok($stats->execute,'execute CoverageStats command '. $stats->command_name);
is(scalar($stats->alignment_summaries),scalar(@wingspans),'found the correcnt number of alignment summaries');
is(scalar($stats->stats_summaries),(scalar(@wingspans)),'found the correct number of stats summaries');
is(scalar($stats->stats_files),(scalar(@wingspans)),'found the correct number of stats files');
