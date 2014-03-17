#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);

my $class = 'Genome::Model::Tools::EpitopePrediction::FrameShift';
my $TEST_DATA_VERSION= 1;
use_ok($class);

my $test_dir = Genome::Utility::Test->data_dir_ok($class, $TEST_DATA_VERSION);
my $bed_file = File::Spec->join($test_dir, "proteome.bed");
my $vcf_file = File::Spec->join($test_dir, "TCGA-A2-A0CM.vcf");
my $database_dir = File::Spec->join($test_dir, "database_dir");
my $expected_output = File::Spec->join($test_dir, "proteome-indel.fasta");
my $expected_output_mod = File::Spec->join($test_dir, "proteome-indel-mod.fasta");

my $output_dir = Genome::Sys->create_temp_directory;

my $cmd = $class->create(
    bed_file => $bed_file,
    vcf_file => $vcf_file,
    database_directory => $database_dir,
    output_directory => $output_dir,
);

ok($cmd->execute, "Command executed");

my $output_file = "$output_dir/proteome-indel.fasta";
my $output_file_mod = "$output_dir/proteome-indel-mod.fasta";
compare_ok($output_file, $expected_output, "Output file is as expected");
compare_ok($output_file_mod, $expected_output_mod, "Output file mod is as expected");

done_testing();
