#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Genome::Model::Tools::Picard::CollectMultipleMetrics;
use Test::More;

my $base_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Picard-CollectMultipleMetrics';
my $bam = $base_dir . '/test.bam';
my $ref = $base_dir . '/test.ref.fa';

my $picard_version = Genome::Model::Tools::Picard->default_picard_version;

my $tmp_dir     = Genome::Sys->create_temp_directory();
my $output_base = $tmp_dir . '/test';

my $cmd = Genome::Model::Tools::Picard::CollectMultipleMetrics->create(
    input_file         => $bam,
    output_basename    => $output_base,
    reference_sequence => $ref,
    program_list       => 'CollectAlignmentSummaryMetrics,CollectInsertSizeMetrics',
    use_version        => $picard_version,
);

ok($cmd, "Command created ok");
ok($cmd->execute, "Command executed ok");

my $ori_align_metrics_file = $base_dir . '/test.alignment_summary_metrics';
my $ori_is_metrics_file    = $base_dir . '/test.insert_size_metrics';

my $align_metrics_file = $output_base . '.alignment_summary_metrics';
my $is_metrics_file    = $output_base . '.insert_size_metrics';

ok(-s $align_metrics_file, 'alignment_summary_metrics generated ok');
ok(-s $is_metrics_file, 'insert_size_metrics generated ok');

my $ori_align_metrics_str = qx(grep -vP '^#' $ori_align_metrics_file);
my $ori_is_metrics_str    = qx(grep -vP '^#' $ori_is_metrics_file);

my $test_align_metrics_str= qx(grep -vP '^#' $align_metrics_file);
my $test_is_metrics_str   = qx(grep -vP '^#' $is_metrics_file);

is($test_align_metrics_str, $ori_align_metrics_str, 'alignment_summary_metrics matches expected');
is($test_is_metrics_str, $ori_is_metrics_str, 'insert_size_metrics matches expected');

done_testing();
