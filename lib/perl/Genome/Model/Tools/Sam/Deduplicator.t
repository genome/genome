#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 2;
use File::Path;


BEGIN
{
    use_ok ('Genome::Model::Tools::Sam::Deduplicator');
}
my ($sam_file1) = ('/gscmnt/sata409/research/mmitreva/edemello/paired/61BNN/6/s_6_1_sequence.txt.PAIRED_REMOVED_1.sam');

my ($deduplicated_file1) = ('/gscmnt/sata409/research/mmitreva/edemello/paired/61BNN/6/s_6_1_sequence.txt.DEDUP_1.sam');

my ($deduplicator1) = (Genome::Model::Tools::Sam::Deduplicator->create(sam_file=>$sam_file1, deduplicated_file=>$deduplicated_file1));

isa_ok($deduplicator1, 'Genome::Model::Tools::Sam::Deduplicator');
