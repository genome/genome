#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 2;
use File::Path;


BEGIN
{
    use_ok ('Genome::Model::Tools::Fasta::Deduplicator');
}

my $fasta_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fasta-Deduplicator/rand60k_short.fna';
my $deduplicated_file = Genome::Sys->create_temp_file_path("dedup.fna");
my $deduplicator = Genome::Model::Tools::Fasta::Deduplicator->create(fasta_file => $fasta_file,
                                                                     deduplicated_file => $deduplicated_file);

my $out = $deduplicator->execute;
ok ($out, "deduplicator runs ok");
