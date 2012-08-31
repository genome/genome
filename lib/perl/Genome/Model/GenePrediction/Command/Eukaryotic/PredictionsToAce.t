#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 13;
use File::Temp 'tempdir';

use_ok('Genome::Model::GenePrediction::Command');
use_ok('Genome::Model::GenePrediction::Command::Eukaryotic');
use_ok('Genome::Model::GenePrediction::Command::Eukaryotic::PredictionsToAce');

my $test_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-GenePredictor/';
ok(-d $test_data_dir, "test data dir exists at $test_data_dir");

my $seq_file = $test_data_dir . '/Contig0a.masked.fasta';
ok(-e $seq_file, "test sequence file exists at $seq_file");

my $expected_ace_output = $test_data_dir . '/example_output.ace';
ok(-e $expected_ace_output, "expected ace output exists at $expected_ace_output");

my $test_output_base_dir = "$ENV{GENOME_TEST_TEMP}/";
ok(-d $test_output_base_dir, "base output dir exists at $test_output_base_dir");

my $temp_dir = tempdir(
    TEMPLATE => 'egap_predictions_to_ace_XXXXXX',
    DIR => $test_output_base_dir,
    CLEANUP => 1,
    UNLINK => 1,
);
ok(-d $temp_dir, "temp output directory exists at $temp_dir");

my $temp_ace_file_fh = File::Temp->new(
    TEMPLATE => 'output_XXXXXX',
    SUFFIX => 'ace',
    DIR => $temp_dir,
    UNLINK => 1,
    CLEANUP => 1,
);
my $temp_ace_file = $temp_ace_file_fh->filename;

my $obj = Genome::Model::GenePrediction::Command::Eukaryotic::PredictionsToAce->create(
    ace_file => $temp_ace_file,
    prediction_directory => $test_data_dir,
    sequence_file => $seq_file,
);
ok($obj, 'successfully created command object');
isa_ok($obj, 'Genome::Model::GenePrediction::Command::Eukaryotic::PredictionsToAce');

ok($obj->execute, 'executed command object');
ok(-s $temp_ace_file, "output written to temp ace file $temp_ace_file");

my $diff_output = Genome::Sys->diff_file_vs_file($temp_ace_file, $expected_ace_output);
ok(!$diff_output, "generated output $temp_ace_file matches expected output $expected_ace_output");
