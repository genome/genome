#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

use above 'Genome';

if ($] < 5.010) {
    plan skip_all => "this test is only runnable on perl 5.10+"
}
plan tests => 6;

map { print STDERR $_ . " => " . $ENV{$_} . "\n" } keys %ENV;

use_ok('Genome::Model::Tools::BioSamtools::CoverageStatsV2');

my $tmp_dir = File::Temp::tempdir('BioSamtools-CoverageStatsV2-'.Genome::Sys->username.'-XXXX',DIR => "$ENV{GENOME_TEST_TEMP}",CLEANUP => 0);

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BioSamtools/RefCov';

my $bam_file = $data_dir .'/test_calmd.bam';
my $regions_file = $data_dir .'/test_regions.bed';
my $genome_file = $data_dir .'/test.genome';

my @wingspans = qw/0 200/;
my @minimum_depths = qw/1 10 20/;
my $stats = Genome::Model::Tools::BioSamtools::CoverageStatsV2->create(
    output_directory => $tmp_dir,
    bam_file => $bam_file,
    bed_file => $regions_file,
    wingspan_values => join(',',@wingspans),
    minimum_depths => join(',',@minimum_depths),
    minimum_base_quality => 20,
    minimum_mapping_quality => 1,
    genome_file => $genome_file,
);
isa_ok($stats,'Genome::Model::Tools::BioSamtools::CoverageStatsV2');
ok($stats->execute,'execute CoverageStats command '. $stats->command_name);
is(scalar($stats->alignment_summaries),scalar(@wingspans),'found the correcnt number of alignment summaries');
is(scalar($stats->stats_summaries),(scalar(@wingspans)),'found the correct number of stats summaries');
is(scalar($stats->stats_files),(scalar(@wingspans)),'found the correct number of stats files');

exit;
