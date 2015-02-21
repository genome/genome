#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Utility::Test 'compare_ok';
use Test::More;
use File::Spec qw();

my $pkg = 'Genome::Model::Tools::Picard::CollectGcBiasMetrics';
use_ok($pkg);

my $data_dir = File::Spec->catfile($ENV{GENOME_TEST_INPUTS}, 'Genome-Model-Tools-Picard-CollectGcBiasMetrics');
my $input_file = File::Spec->catfile($data_dir, "sorted.bam");
my $ref_file = File::Spec->catfile($data_dir, "small.fa");
my $expected_file = File::Spec->catfile($data_dir, "expected.txt");

my $output_file = Genome::Sys->create_temp_file_path;

my $cmd = $pkg->create(
    input_file => $input_file,
    refseq_file => $ref_file,
    output_file => $output_file,
    );

ok($cmd, "created command");
ok($cmd->execute, "executed command");

compare_ok($expected_file, $output_file, filters => ['^#.*']);

done_testing();
