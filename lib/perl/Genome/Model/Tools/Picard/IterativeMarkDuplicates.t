#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More tests => 4;
use File::Compare;

use_ok('Genome::Model::Tools::Picard::IterativeMarkDuplicates');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Picard-IterativeMarkDuplicates';
my $test_bam = $data_dir .'/test.bam';
my $expected_metrics = $data_dir .'/expected_metrics.tsv';

my $tmp_dir = File::Temp::tempdir('Picard-IterativeMarkDuplicates-'.Genome::Sys->username.'-XXXX',DIR => "$ENV{GENOME_TEST_TEMP}",CLEANUP => 1);

my $iterative_mrkdup = Genome::Model::Tools::Picard::IterativeMarkDuplicates->create(
    input_file => $test_bam,
    output_directory => $tmp_dir,
    step_size => 1000,
);
isa_ok($iterative_mrkdup,'Genome::Model::Tools::Picard::IterativeMarkDuplicates');
ok($iterative_mrkdup->execute,'execute command '. $iterative_mrkdup->command_name);
ok(!compare($expected_metrics,$iterative_mrkdup->output_file),'Expected BAM '. $expected_metrics
       .' is identical to IterativeMarkDuplicates metrics '. $iterative_mrkdup->output_file);
#ok(unlink($iterative_mrkdup->output_file),'Removed temporary output file : '. $iterative_mrkdup->output_file);

exit;
