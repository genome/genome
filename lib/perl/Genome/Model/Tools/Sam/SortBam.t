#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Genome::Model::Tools::Sam::SortBam;
use Test::More;
#tests => 1;

if (`uname -a` =~ /x86_64/){
    plan tests => 3;
} else{
    plan skip_all => 'Must run on a 64 bit machine';
}

my $input = $ENV{GENOME_TEST_INPUTS} . '/Genome-Tools-Sam-SortBam/normal.tiny.bam';

# step 1: test 1 file case

my $out_file = File::Temp->new(SUFFIX => ".bam" );

my $cmd_1 = Genome::Model::Tools::Sam::SortBam->create(file_name=>$input,
                                                     output_file=>$out_file->filename);

ok($cmd_1, "created command");
ok($cmd_1->execute, "executed");
ok(-s $out_file->filename, "output file is nonzero");


print "done\n";

