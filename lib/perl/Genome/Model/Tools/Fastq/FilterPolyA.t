#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";  # forces a 'use lib' when run directly from the cmdline
use Test::More tests => 3;
use FindBin qw($Bin);
use File::Temp;
 
my $in = $Bin . '/' . 'FilterPolyA.t.d/in.fastq';

my $dir = File::Temp::tempdir(CLEANUP => 1);
$dir or die "Failed to create temp directory!";
my $out = "$dir/out.fastq";

my $result = Genome::Model::Tools::Fastq::FilterPolyA->execute(input1 => $in, output1 => $out);
ok($result, "got result");

ok(-e $out, "ouput file exists");

my $expected_out = "$Bin/FilterPolyA.t.d/expected-out.fastq";

my $differences = `diff $out $expected_out | wc -l`;
chomp $differences;

is($differences,0,"no difference between actual and expected output");

