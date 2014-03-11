#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);

my $TEST_DATA_VERSION = 1;
my $class = 'Genome::Model::Tools::EpitopePrediction::GenerateVariantSeq';
use_ok($class);

my $test_dir = Genome::Utility::Test->data_dir_ok($class, $TEST_DATA_VERSION);
my $input_file = File::Spec->join($test_dir, "input.tsv");
for my $length (qw(17 21 31)) {
    test_for_length($length, $test_dir, $input_file);
}

sub test_for_length {
    my ($length, $test_dir, $input_file) = @_;

    my $expected_output = File::Spec->join($test_dir, "output_" . $length ."mer");
    my $output_file = Genome::Sys->create_temp_file_path;

    my $cmd = $class->create(
        input_file => $input_file,
        output_file => $output_file,
        length => $length,
    );
    ok($cmd, "Created a command for length $length");

    ok($cmd->execute, "Command executed for length $length");

    compare_ok($output_file, $expected_output, "Output file is as expected for length $length");
}

done_testing();
