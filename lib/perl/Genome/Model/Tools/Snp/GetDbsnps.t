#!/usr/bin/env genome-perl

use strict;
use warnings;
use above "Genome";
use File::Temp;
use Test::More tests => 5;

my $chromosome = "10";
my $start = 126008345;
my $stop = 126010576;
my $organism = "human"; # or mouse

my $out_fh = File::Temp->new(
    TEMPLATE => 'Genome-Model-Tools-Snp-GetDbsnps-XXXXXX',
    DIR => "$ENV{GENOME_TEST_TEMP}/",
    CLEANUP => 1,
    UNLINK => 1,
);
my $out = $out_fh->filename;
$out_fh->close;

my $dbsnpout = $out . ".dbsnp.gff";

my $expected_output = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Snp-GetDbsnps/GetDbsnps.expectedout.dbsnp.gff";

ok(-e $expected_output, "expected output found at $expected_output");

my $dbsnps = Genome::Model::Tools::Snp::GetDbsnps->create(chromosome=>$chromosome,start=>$start,stop=>$stop,organism=>$organism,gff=>1,out=>$out);
ok($dbsnps, "command object created");

ok($dbsnps->execute(), "executed command");
ok(-e $dbsnpout, "dbsnp output file exists");

my $diff = `diff $dbsnpout $expected_output`;
ok($diff eq '', "output as expected");

