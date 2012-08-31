#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}
else {
    plan tests => 32;
}

use above 'Genome';

use_ok('Genome::Model::Tools::Dbsnp::Import');

my $sub_dir = "/Genome-Model-Tools-Dbsnp-Import";
my $test_dir = $ENV{GENOME_TEST_INPUTS}.$sub_dir;
ok (-d $test_dir, "test directory $test_dir is present");
my $test_output = "$test_dir/v2/output.bed";
ok(-e $test_output, "test output $test_output exists");
my $test_url = $ENV{GENOME_TEST_URL}.$sub_dir;

my @chromosomes = qw( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT Un);
for my $chromosome (@chromosomes){
    my $flatfile = "$test_dir/ds_flat_ch$chromosome.flat";
    ok(-e $flatfile, "test file $flatfile exists");
}

my $command_output = Genome::Sys->create_temp_file_path();

my $cmd = Genome::Model::Tools::Dbsnp::Import->create(
    output_file => $command_output,
    flat_file_url => $test_url,
    filename_pattern => "ds_flat_chX.flat",
);
ok ($cmd, 'created the importer');
ok($cmd->execute, 'importer ran successfully');

#This is a hack, since I can't figure out how to make Joinx sort past chromosome, start, and stop
my $sorted_command_output = Genome::Sys->create_temp_file_path();
system("sort $command_output > $sorted_command_output");
my $sorted_test_output = Genome::Sys->create_temp_file_path();
system("sort $test_output > $sorted_test_output");

my $diff = Genome::Sys->diff_file_vs_file($sorted_test_output, $sorted_command_output);
ok(!$diff, 'returned file matches expected file')
    or diag("diff:\n" . $diff);
