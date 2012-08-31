#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 1;
use File::Path;

BEGIN {
        use_ok ('Genome::Model::Tools::Fastq::TrimPrimer');
}

my $trim = Genome::Model::Tools::Fastq::TrimPrimer->create(fastq_file=>'/gscmnt/sata413/research/kwylie/2010_03_08_virome_pipeline_test_runs/Illumina_4/s_4_1_sequence.txt', 
                                                           output=>'/gscmnt/sata413/research/kwylie/2010_03_08_virome_pipeline_test_runs/Illumina_4/FOO.txt');

exit;

