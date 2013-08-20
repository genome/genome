#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Genome::Utility::Test qw(compare_ok);
use Test::More tests => 9;

my $class = 'Genome::Model::Tools::Varscan::ConsensusVcfMatch';

use_ok($class);

my $test_dir = Genome::Utility::Test->data_dir($class);
$test_dir .= "/v1";

my $input_tsv = $test_dir . "/input.tsv";
my $input_cns = $test_dir . "/input.cns";
my $expected_output_tsv = $test_dir . "/output.tsv";
my $expected_output_cns = $test_dir . "/output.cns";

ok(-e $input_tsv, "test file $input_tsv exists");
ok(-e $input_cns, "test file $input_cns exists");
ok(-e $expected_output_tsv, "test file $expected_output_tsv exisits");
ok(-e $expected_output_cns, "test file $expected_output_cns exisits");

my $test_output_tsv = Genome::Sys->create_temp_file_path();
my $test_output_cns = Genome::Sys->create_temp_file_path();

my $cmd = Genome::Model::Tools::Varscan::ConsensusVcfMatch->create(
    indel_vcf => $input_tsv,
    sample_cns_list_file => $input_cns,
    output_file => $test_output_tsv,
    output_cns_file => $test_output_cns,
);

ok($cmd, "Successfully created command");
ok($cmd->execute, "Successfully ran command");

my $tsv_diff = `diff -q $expected_output_tsv $test_output_tsv`;
ok($tsv_diff eq '' , "output file matches expected output file");

my $cns_diff = `diff -q $expected_output_cns $test_output_cns`;
ok($cns_diff eq '' , "output cns file matches expected output cns file");
