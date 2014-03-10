#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test;

my $class = 'Genome::Model::SmallRna::Command::StatsGenerator';
use_ok($class);

my $data_dir    = Genome::Utility::Test->data_dir_ok($class);
my $test_dir    = Genome::Sys->create_temp_directory;
my $exp_out_dir = $data_dir . '/v1';

my $bam_file            = $data_dir.'/test.bam';
my $stats_file          = $data_dir.'/coverage_stats.tsv';
my $normalized_flagstat = $data_dir.'/17_70.bam.flagstat';

my @out_file_names = qw(alignment_stats.tsv top_sorted_clusters.bed subclusters.bed sub_inter_file);
my @test_files = map{$test_dir.'/'.$_}@out_file_names;
my @exp_files  = map{$exp_out_dir.'/'.$_}@out_file_names;

my $cmd = $class->execute(
    bam_file             => $bam_file,
    coverage_stats_file  => $stats_file,
    flagstat_17_70_file  => $normalized_flagstat,
    output_stats_file    => $test_files[0],
    output_clusters_file => $test_files[1],
    output_subclusters_file          => $test_files[2],
    output_subcluster_intersect_file => $test_files[3],
);

#skip sub_inter_file since no output from the test data
for my $i (0..2) {
    Genome::Utility::Test::compare_ok($test_files[$i], $exp_files[$i], 'output file '.$out_file_names[$i].' created ok');
}

done_testing;

