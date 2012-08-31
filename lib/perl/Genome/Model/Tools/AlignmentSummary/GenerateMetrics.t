#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Compare;
use Test::More tests => 4;

use_ok('Genome::Model::Tools::AlignmentSummary::GenerateMetrics');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-AlignmentSummary-GenerateMetrics';

my $test_bed = $data_dir .'/22.bed';
my $test_bam = $data_dir .'/22.bam';
my $expected_as = $data_dir .'/expected_as.yaml';

my $tmp_dir = File::Temp::tempdir('AlignmentSummary-GenerateMetrics-'.Genome::Sys->username.'-XXXX',DIR => "$ENV{GENOME_TEST_TEMP}",CLEANUP => 1);
my $output_as = $tmp_dir .'/as.yaml';

my $as_cmd = Genome::Model::Tools::AlignmentSummary::GenerateMetrics->create(
    use_version => '1.2.6',
    alignment_file_path => $test_bam,
    roi_file_path => $test_bed,
    output_file => $output_as,
);
isa_ok($as_cmd,'Genome::Model::Tools::AlignmentSummary::GenerateMetrics');
ok($as_cmd->execute,'execute command '. $as_cmd->command_name);

ok( (compare($as_cmd->output_file,$expected_as) == 0),'The alignment summary YAML file is the same!');

exit;
