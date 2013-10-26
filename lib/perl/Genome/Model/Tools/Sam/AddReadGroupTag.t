#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Genome::Model::Tools::Sam::Merge;
use Test::More tests => 5;

my $input = $ENV{GENOME_TEST_INPUTS} . "/Genome-InstrumentData-Alignment/test.sam";

# step 1: test 1 
my $tmp_dir = File::Temp->newdir( "AddReadGroupTag_XXXXX",
                                  TMPDIR => 1,
                                  CLEANUP => 1 );


my $output_file = File::Temp->new(SUFFIX => ".sam", DIR => $tmp_dir);

#uncomment to inspect output 
$output_file->unlink_on_destroy(1);

my $cmd_1 = Genome::Model::Tools::Sam::AddReadGroupTag->create(input_file=>$input,
                                                              output_file=>$output_file->filename,
                                                              read_group_tag=>"123456",
                                                            );

ok($cmd_1, "created command");
ok($cmd_1->execute, "executed");
ok(-s $output_file->filename, "output file is nonzero");

open(FILE,$output_file) or die("Unable to open ".$output_file->filename. " in AddReadGroupTag.t");
my @result_file = <FILE>;
close(FILE);

my @rg_result = grep (/RG:Z:123456/,@result_file);
my @pg_result = grep (/PG:Z:123456/,@result_file);

ok( scalar(@rg_result) eq 49, "found 49 RG tags.");
ok( scalar(@pg_result) eq 49, "found 49 PG tags.");
