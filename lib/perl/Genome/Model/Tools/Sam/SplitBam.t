#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Genome::Model::Tools::Sam::SplitBam;
use Test::More;
#tests => 1;

if (`uname -a` =~ /x86_64/){
    plan tests => 3;
} else{
    plan skip_all => 'Must run on a 64 bit machine';
}

my $bam_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Tools-Sam-SortBam/normal.tiny.bam';
my $size = 1000;
my $expected_number_bam_files = 6;
my $tmp_out_dir = File::Temp::tempdir(CLEANUP => 1);

my $cmd_1 = Genome::Model::Tools::Sam::SplitBam->create(
    bam_file => $bam_file,
    output_directory => $tmp_out_dir,
    size => $size,
);

ok($cmd_1, 'created command');
ok($cmd_1->execute, 'executed');
my $sub_bam_files = $cmd_1->sub_bam_files;
is(scalar(@{$sub_bam_files}),$expected_number_bam_files,'Found the expected number of sub bam files');
