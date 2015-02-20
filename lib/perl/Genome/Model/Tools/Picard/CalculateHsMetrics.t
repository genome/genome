#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);
use File::Spec qw();

my $pkg = 'Genome::Model::Tools::Picard::CalculateHsMetrics';

use_ok($pkg);

my $data_dir = File::Spec->catfile($ENV{GENOME_TEST_INPUTS}, 'Genome-Model-Tools-Picard-CalculateHsMetrics');
my $bam_file = File::Spec->catfile($data_dir, 'coordsort.bam');
my $expected_file = File::Spec->catfile($data_dir, 'expected.txt');
my $targets_file = File::Spec->catfile($data_dir, 'targets.txt');
my $output_file = Genome::Sys->create_temp_file_path;

my $obj = $pkg->create(
    input_file => $bam_file,
    output_file => $output_file,
    target_intervals => $targets_file,
    bait_intervals => $targets_file,
    );
ok($obj, 'created command');
ok($obj->execute, 'executed command');

compare_ok($expected_file, $output_file, diag => 1, filters => ['^#.*']);
done_testing();
