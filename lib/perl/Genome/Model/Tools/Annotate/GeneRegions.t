#!/usr/bin/env genome-perl

use strict;
use warnings;
use above "Genome";
use File::Temp;
use Test::More tests => 6;

my $list = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Annotate-GeneRegions/coordinates_list";
ok(-f $list, "coordinates list file found at $list");

my $output_fh = File::Temp->new(
    TEMPLATE => 'Genome-Model-Tools-Annotate-GeneRegions-XXXXXX',
    DIR => "$ENV{GENOME_TEST_TEMP}/",
    CLEANUP => 1,
    UNLINK => 1,
);
my $output = $output_fh->filename;
$output_fh->close;

my $regions = Genome::Model::Tools::Annotate::GeneRegions->create(list=>$list, mode => 'coordinates', output => $output);
ok($regions, "GeneRegions object successfully created");
ok($regions->execute(), "successfully executed command");
ok(-f $output, "object file exists");
my $expected_output = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Annotate-GeneRegions/expected_output";
ok(-f $expected_output, "expected output file exists");

my $diff = `diff $output $expected_output`;
ok($diff eq '', "output as expected");


