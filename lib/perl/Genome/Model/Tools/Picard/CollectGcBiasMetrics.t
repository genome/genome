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

my $tmpdir = Genome::Sys->create_temp_directory;

my $output_file = File::Spec->catfile($tmpdir, 'gc_bias.txt');
my $summary_output_file = File::Spec->catfile($tmpdir, 'gc_bias_summary.txt');
my $chart_output_file = File::Spec->catfile($tmpdir, 'gc_bias.pdf');

my $cmd = $pkg->create(
    input_file => $input_file,
    refseq_file => $ref_file,
    output_file => $output_file,
    summary_output => $summary_output_file,
    chart_output => $chart_output_file,
    );

ok($cmd, "created command");
ok($cmd->execute, "executed command");

compare_ok($expected_file, $output_file, filters => ['^#.*']);

done_testing();
