#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 2;
use File::Path;


BEGIN
{
    use_ok ('Genome::Model::Tools::Sam::RemovePairedEnds');
}

my $sam_file_1 =    '/gscmnt/sata409/research/mmitreva/edemello/paired/61BKE/2/s_2_1_sequence.txt.UNALIGNED_1.sam';
my $sam_file_2 =    '/gscmnt/sata409/research/mmitreva/edemello/paired/61BKE/2/s_2_2_sequence.txt.UNALIGNED_2.sam';
my $paired1 =       '/gscmnt/sata409/research/mmitreva/edemello/paired/61BKE/2/s_2_1_sequence.txt.PAIRED.sam';
my $paired2 =       '/gscmnt/sata409/research/mmitreva/edemello/paired/61BKE/2/s_2_2_sequence.txt.PAIRED.sam';

my $per = Genome::Model::Tools::Sam::RemovePairedEnds->create(  sam1                        => $sam_file_1,
                                                                sam2                        => $sam_file_2,
                                                                paired_end_removed_file1    => $paired1,
                                                                paired_end_removed_file2    => $paired2,);

isa_ok($per, 'Genome::Model::Tools::Sam::RemovePairedEnds');
