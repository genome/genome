#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Utility::Test 'compare_ok';
use Test::More;
use Test::Exception;
use File::Spec qw();

my $pkg = 'Genome::Model::Tools::Picard::BamIndexStats';
use_ok($pkg);

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Picard-BamIndexStats';

my $input_file = File::Spec->catfile($data_dir, "coordsort.bam");
my $expected_file = File::Spec->catfile($data_dir, "expected.txt");
my $output_file = Genome::Sys->create_temp_file_path;

my $cmd = $pkg->create(
    input_file => "$data_dir/coordsort.bam",
    output_file => $output_file
    );

ok($cmd, "created command");
ok($cmd->execute, "executed command");

compare_ok($expected_file, $output_file, "output file is correct");

my $bad_cmd = $pkg->create(
    input_file => "$data_dir/coordsort.bam",
    output_file => $output_file,
    use_version => "1.22",
    );

ok($bad_cmd, "created command with invalid version");
dies_ok(sub { $bad_cmd->execute }, "command with invalid version fails to execute");

done_testing();
