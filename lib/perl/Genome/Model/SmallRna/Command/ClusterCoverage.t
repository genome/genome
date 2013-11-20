#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test;

my $class = 'Genome::Model::SmallRna::Command::ClusterCoverage';
use_ok($class);

plan skip_all => 'This test takes more than 10 mins to run. Skip for now';

my $data_dir    = Genome::Utility::Test->data_dir_ok($class);
my $test_dir    = Genome::Sys->create_temp_directory;
my $exp_out_dir = $data_dir . '/v1';

for my $num (1, 2) {
    run_test($data_dir, $test_dir, $num);
}

done_testing;


sub run_test {
    my ($data_dir, $test_dir, $num) = @_;

    my $bam_file   = $data_dir.'/test_'.$num.'.bam';
    my $stats_file = $test_dir.'/stats_file'.$num;
    my $bed_file   = $test_dir.'/bed_file'.$num;

    my $exp_stats_file = $exp_out_dir . '/stats_file'.$num;
    my $exp_bed_file   = $exp_out_dir . '/bed_file'.$num;

    my $cmd = $class->execute(
        bam_file       => $bam_file,
        zenith_depth   => 5,
        minimum_depth  => 1,
        stats_file     => $stats_file,
        bed_file       => $bed_file,
    );

    Genome::Utility::Test::compare_ok($stats_file, $exp_stats_file, "Output stats file $num was created as expected");
    Genome::Utility::Test::compare_ok($bed_file, $exp_bed_file, "Output bed file $num was created as expected");
}

