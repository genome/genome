#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Genome::Model::Tools::Picard::MergeSamFiles;
use Test::More tests => 6;

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Tools-Sam-Merge';
my $input_normal = $dir. '/normal.tiny.bam';
my $input_tumor  = $dir. '/tumor.tiny.bam';

# step 1: test 1 file case

my $out_1_file = File::Temp->new(SUFFIX => ".bam" );

my $cmd_1 = Genome::Model::Tools::Picard::MergeSamFiles->create(
    input_files => [$input_normal],
    output_file    => $out_1_file->filename,
);

ok($cmd_1, "created command");
ok($cmd_1->execute, "executed");
ok(-s $out_1_file->filename, "output file is nonzero");

# step 1: test >1 input file case

my $out_2_file = File::Temp->new(SUFFIX => ".bam" );

my $cmd_2 = Genome::Model::Tools::Picard::MergeSamFiles->create(
    input_files => [$input_normal, $input_tumor],
    output_file    => $out_2_file->filename,
);

ok($cmd_2, "created command");
ok($cmd_2->execute, "executed");
ok(-s $out_2_file->filename, "output file is nonzero");


