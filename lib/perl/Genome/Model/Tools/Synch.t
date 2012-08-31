#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 2;
use File::Path;


BEGIN
{
    use_ok ('Genome::Model::Tools::Synch');
}

my $fastq_file   = '/gscmnt/sata156/research/mmitreva/bacteria/TwoVisit_TwoPeople/data/61BKE/s_2_1_sequence.txt';
my $fastq_file2  = '/gscmnt/sata156/research/mmitreva/bacteria/TwoVisit_TwoPeople/data/61BKE/s_2_2_sequence.txt';
my $name1        = '/gscmnt/sata409/research/mmitreva/edemello/paired/61BKE/2/s_2_1_sequence.txt.RESURRECTED_1.sam';
my $name2        = '/gscmnt/sata409/research/mmitreva/edemello/paired/61BKE/2/s_2_2_sequence.txt.RESURRECTED_2.sam';
my $output       = 'SYNCH';
my $prefix1      = '/gscmnt/sata409/research/mmitreva/edemello/paired/61BKE/2/s_2_1_sequence.txt.';
my $prefix2      = '/gscmnt/sata409/research/mmitreva/edemello/paired/61BKE/2/s_2_2_sequence.txt.';


my $synch = Genome::Model::Tools::Synch->create(  name1          => $name1,
                                                  name2          => $name2,
                                                  fastq1         => $fastq_file,
                                                  fastq2         => $fastq_file2,
                                                  output         => $output,
                                                  prefix1        => $prefix1,
                                                  prefix2        => $prefix2,
);

isa_ok($synch, 'Genome::Model::Tools::Synch');
