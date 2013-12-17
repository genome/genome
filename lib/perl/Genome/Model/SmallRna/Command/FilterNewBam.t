#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test;

plan skip_all => 'This test takes 7 mins to run. Skip for now';

my $class = 'Genome::Model::SmallRna::Command::FilterNewBam';
use_ok($class);

my $data_dir  = Genome::Utility::Test->data_dir_ok($class);
my $test_file = Genome::Sys->create_temp_file_path;

my $bam_file  = $data_dir.'/test.bam';

my $exp_out_dir  = $data_dir . '/v1';
my $exp_out_file = $exp_out_dir.'/61_70.bam';

my $cmd = $class->execute(
    bam_file          => $bam_file,
    filtered_bam_file => $test_file,
    read_size_bin     => '61_70',
    xa_tag            => 1,
);

Genome::Utility::Test::compare_ok($test_file, $exp_out_file, 'output file 61_70.bam created ok');

done_testing;

