#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);

my $class = 'Genome::Model::Tools::EpitopePrediction::RunNetmhcNew';
my $TEST_DATA_VERSION= 1;
use_ok($class);

my $test_dir = Genome::Utility::Test->data_dir_ok($class, $TEST_DATA_VERSION);
my $input_fasta = File::Spec->join($test_dir, "snvs_21mer.fa");
my $expected_output = File::Spec->join($test_dir, "expected.xls");
my $output_file = Genome::Sys->create_temp_file_path;
my $stdout_file = Genome::Sys->create_temp_file_path;

my $cmd = $class->create(
    allele => 'HLA-A02:01',
    fasta_file => $input_fasta,
    output_file => $output_file,
    epitope_length => 9,
    stdout_file => $stdout_file,
);
ok($cmd, "Created command");

ok($cmd->execute, "Command executed");

compare_ok($output_file, $expected_output, filters => ['^NetMHC version.*']);

done_testing();
