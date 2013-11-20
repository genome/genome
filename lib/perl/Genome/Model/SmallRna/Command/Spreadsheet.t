#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test;

my $class = 'Genome::Model::SmallRna::Command::Spreadsheet';
use_ok($class);

my $data_dir    = Genome::Utility::Test->data_dir_ok($class);
my $test_file   = Genome::Sys->create_temp_file_path;

my $annot_inter_file = $data_dir.'/annotation_intersect.tsv';
my $stat_file        = $data_dir.'/alignment_stats.tsv';

my $exp_out_dir  = $data_dir . '/v1';
my $exp_out_file = $exp_out_dir.'/Final_spreadsheet.tsv';

my $cmd = $class->execute(
    input_stats_file     => $stat_file,
    input_intersect_file => $annot_inter_file,
    output_spreadsheet   => $test_file,
    input_cluster_number => 5000,
);

Genome::Utility::Test::compare_ok($test_file, $exp_out_file, 'output file Final_spreadsheet.tsv created ok');

done_testing;

