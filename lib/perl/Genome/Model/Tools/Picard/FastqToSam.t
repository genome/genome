#!/usr/bin/env genome-perl

use strict;
use warnings FATAL => 'all';

use above 'Genome';
use Genome::Model::Tools::Picard::FastqToSam;
use Test::More tests => 3;
use File::Temp;
use Path::Class qw(dir file);

# data here is first 100 lines from lane 1 of
# /gscmnt/sata604/hiseq2000/100218_P21_0393_AFC20GF1/Data/Intensities/Basecalls/GERALD_30-03-2010_lims
my $dir = dir(
    $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Picard-FastqToSam');

# step 1: test 1 file case

my $out_1_file = File::Temp->new( SUFFIX => '.bam' );
my $cmd_1 = Genome::Model::Tools::Picard::FastqToSam->create(
    fastq           => $dir->file('s_1_1_sequence.txt') . '',
    fastq2          => $dir->file('s_1_2_sequence.txt') . '',
    output          => $out_1_file->filename,
    platform_unit   => 'AFC20GF1.1',
    quality_format  => 'Illumina',
    platform        => 'Illumina',
    sample_name     => 'H_KU-6888-D59687',
    read_group_name => 2854142190,
    sort_order      => 'queryname',
    library_name    => 'foo',
);
isa_ok( $cmd_1, 'Genome::Model::Tools::Picard::FastqToSam' );
ok( $cmd_1->execute,          'execute' );
ok( -s $out_1_file->filename, 'output file is non-zero' );
