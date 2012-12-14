#!/usr/bin/env genome-perl

use strict;
use warnings;
use above "Genome";
use File::Temp;
use Test::More tests => 4;

my $chromosome = "1";	
my $start = 1;
my $stop = 10000;
my $organism = "human";
my $version = "54_36p_v2";

my $output_fh = File::Temp->new(
    TEMPLATE => 'Genome-Model-Tools-Annotate-TranscriptRegions-XXXXXX',
    DIR => "$ENV{GENOME_TEST_TEMP}",
    CLEANUP => 1,
    UNLINK => 1,
);
my $output = $output_fh->filename;
$output_fh->close;

my $expected_output = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Annotate-TranscriptRegions/TranscriptRegions.t.expected.output.txt.new";
ok(-e $expected_output, "expected output file exists at $expected_output");

my $regions = Genome::Model::Tools::Annotate::TranscriptRegions->create(chromosome=>$chromosome,start=>$start,stop=>$stop,organism=>$organism,version=>$version,output=>$output);
ok($regions, "TranscriptRegions command object created");

ok($regions->execute(), "successfully executed command");

my $diff = `sdiff -s $output $expected_output`;
ok($diff eq '', "output as expected") or diag($diff);
