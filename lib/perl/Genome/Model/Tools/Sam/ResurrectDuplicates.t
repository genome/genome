#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 2;
use File::Path;


BEGIN
{
    use_ok ('Genome::Model::Tools::Sam::ResurrectDuplicates');
}

my $raw_dedup_sam1 = '/gscmnt/sata409/research/mmitreva/edemello/paired/61BNN/6/s_6_1_sequence.txt.DEDUP_1.sam';
my $raw_dedup_sam2 = '/gscmnt/sata409/research/mmitreva/edemello/paired/61BNN/6/s_6_2_sequence.txt.DEDUP_2.sam';
my $orig_sam1 = '/gscmnt/sata409/research/mmitreva/edemello/paired/61BNN/6/s_6_1_sequence.txt.UNALIGNED_1.sam'; 
my $orig_sam2 = '/gscmnt/sata409/research/mmitreva/edemello/paired/61BNN/6/s_6_2_sequence.txt.UNALIGNED_2.sam'; 
my $output_file1 = '/gscmnt/sata409/research/mmitreva/edemello/paired/61BNN/6/s_6_1_sequence.txt.RESURRECTED_1.sam';  
my $output_file2 = '/gscmnt/sata409/research/mmitreva/edemello/paired/61BNN/6/s_6_2_sequence.txt.RESURRECTED_2.sam';

my $rd = Genome::Model::Tools::Sam::ResurrectDuplicates->create(            
            raw_dedup_sam1 => $raw_dedup_sam1, 
            raw_dedup_sam2 => $raw_dedup_sam2, 
            orig_sam1 => $orig_sam1, 
            orig_sam2 => $orig_sam2, 
            output_file1 => $output_file1,     
            output_file2 => $output_file2,     
);

isa_ok($rd, 'Genome::Model::Tools::Sam::ResurrectDuplicates');
