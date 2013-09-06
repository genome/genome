#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

use above 'Genome';

use_ok('Genome::Model::Tools::Mutect::MergeOutputFiles');
my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Mutect-Merge-Output-Files";

my $expected_file = "$test_dir/expected.v1";
my @input_files = ("$test_dir/output1.v1", "$test_dir/output2.v1");

my $output_file = Genome::Sys->create_temp_file_path;
my $command = Genome::Model::Tools::Mutect::MergeOutputFiles->create(
    mutect_output_files => \@input_files,
    merged_file => $output_file,
);
ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file created");

my $diff = Genome::Sys->diff_file_vs_file($output_file, $expected_file);
    ok(!$diff, 'output matched expected result')
        or diag("diff results:\n" . $diff);

done_testing();

