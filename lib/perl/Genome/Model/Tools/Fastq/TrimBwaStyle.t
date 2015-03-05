#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 5;
use File::Compare;

use above 'Genome';

use_ok('Genome::Model::Tools::Fastq::TrimBwaStyle') or die;

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $base_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fastq-TrimBwaStyle';
my $fastq_file = "$base_dir/test.fastq";

my $trim = Genome::Model::Tools::Fastq::TrimBwaStyle->create(
    fastq_file  => $fastq_file,
    out_file    => $tmp_dir.'/test.trimmed.fastq',
);
isa_ok($trim,'Genome::Model::Tools::Fastq::TrimBwaStyle');

ok($trim->execute,'execute command '. $trim->command_name);

for my $file (qw(test.trimmed.fastq trim.report)) {
    my $output_file = $tmp_dir."/$file";
    my $expect_file = $base_dir."/$file";
    ok(compare($output_file, $expect_file) == 0, "Output: $file is created as expected");
}
