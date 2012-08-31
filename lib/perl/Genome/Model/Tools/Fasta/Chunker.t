#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use File::Compare;
use Test::More tests => 8;


BEGIN {use_ok('Genome::Model::Tools::Fasta::Chunker');}

my ($dir, $tmp_dir, $chunk_size) = ($ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fasta-Chunker/', Genome::Sys->create_temp_directory, 5);

#create
my $chunker = Genome::Model::Tools::Fasta::Chunker->create(
                fasta_file=> $dir . 'static.fna',
                tmp_dir => $tmp_dir,
                chunk_size=> $chunk_size,
);
isa_ok($chunker, 'Genome::Model::Tools::Fasta::Chunker');

#chunk file
ok($chunker->execute, 'chunking file');

#compare to static chunks 
my $file_chunks = $chunker->file_chunks;

for (my ($i,$static_chunk) = (0,undef); $i < $chunk_size; $i++)
{
    $static_chunk = $dir . 'CHUNK_' . $i . '.static';    
    cmp_ok(compare(@$file_chunks[$i],$static_chunk), '==', 0, "$static_chunk matches " . @$file_chunks[$i]);
}

