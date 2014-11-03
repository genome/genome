#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";  # forces a 'use lib' when run directly from the cmdline
use Test::More tests => 4;
use File::Temp;

my $data_dir = __FILE__ . '.d';
my $in = "$data_dir/in.fastq";

my $dir = File::Temp::tempdir(CLEANUP => 1);
$dir or die "Failed to create temp directory!";
my $out = "$dir/out.fastq";

my $cmd = Genome::Model::Tools::Fastq::FilterPolyA->execute(input1 => $in, output1 => $out);
ok($cmd, "got command");
ok($cmd->result, "got result");

ok(-e $out, "ouput file exists");

my $expected_out = "$data_dir/expected-out.fastq";

my $differences = `diff $out $expected_out | wc -l`;
chomp $differences;

is($differences,0,"no difference between actual and expected output");
