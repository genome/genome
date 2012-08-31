#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Model::Tools::Picard::SortSam;
use Test::More tests => 6;

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Tools-Sam-Merge';
my $input_normal = $dir. '/normal.tiny.bam';
my $input_tumor  = $dir. '/tumor.tiny.bam';

# step 1: test 1 file case

my $out_1_file = File::Temp->new(SUFFIX => '.bam' );
my $cmd_1 = Genome::Model::Tools::Picard::SortSam->create(
    input_file => $input_normal,
    output_file    => $out_1_file->filename,
    use_version => 'r107',
);
isa_ok($cmd_1,'Genome::Model::Tools::Picard::SortSam');
my $rv;
eval {
    $rv = $cmd_1->execute
};
ok(!$rv, 'correctly failed to execute using an old version of Picard');
ok(-z $out_1_file->filename, 'output file is zero');

# test 2: use an updated version of picard

my $out_2_file = File::Temp->new(SUFFIX => '.bam' );
my $cmd_2 = Genome::Model::Tools::Picard::SortSam->create(
    input_file => $input_normal,
    output_file    => $out_2_file->filename,
    use_version => '1.17',
);
isa_ok($cmd_2,'Genome::Model::Tools::Picard::SortSam');
ok($cmd_2->execute, 'execute');
ok(-s $out_2_file->filename, 'output file is non-zero');



