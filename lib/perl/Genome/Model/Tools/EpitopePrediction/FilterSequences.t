#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Test::Exception;
use Genome::Utility::Test qw(compare_ok);

my $class = 'Genome::Model::Tools::EpitopePrediction::FilterSequences';
my $TEST_DATA_VERSION= 3;
use_ok($class);

my $test_dir = Genome::Utility::Test->data_dir_ok($class, $TEST_DATA_VERSION);

subtest "filter ok" => sub {
    for my $file (qw(input_unknown_sequences.fasta input_duplicate_sequences.fasta)) {
        my $input_file = File::Spec->join($test_dir, $file);
        my $expected_output = File::Spec->join($test_dir, "expected.fasta");
        my $output_dir = Genome::Sys->create_temp_directory;

        my $cmd = $class->create(
            input_file => $input_file,
            output_directory => $output_dir,
        );

        ok($cmd->execute, sprintf("Command with input file (%s) executed", $file));

        compare_ok($cmd->output_file, $expected_output, "Output file is as expected");
    }
};

subtest "error" => sub {
    my $file = 'input_duplicate_id_different_sequence.fasta';
    my $input_file = File::Spec->join($test_dir, $file);
    my $output_dir = Genome::Sys->create_temp_directory;

    my $cmd = $class->create(
        input_file => $input_file,
        output_directory => $output_dir,
    );

    throws_ok(sub {$cmd->execute}, qr/Found duplicate entries with id \(MT\.Rhbdd3\.p\.Q191L\) but different sequences/, "Command with input file ($file) fails as expected");
};

done_testing();
