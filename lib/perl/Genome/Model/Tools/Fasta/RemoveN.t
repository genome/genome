#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 2;
use File::Path;


BEGIN
{
    use_ok ('Genome::Model::Tools::Fasta::RemoveN');
}

my $fasta_file = '/gscmnt/temp206/info/seqana/species_independant/edemello/short.fna';
my $n_removed_file = '/gscmnt/temp206/info/seqana/species_independant/edemello/illumina/illumina_bwa.N_REMOVED.fasta';

my $n_remover = Genome::Model::Tools::Fasta::RemoveN->create(fasta_file     =>  $fasta_file,
                                                             n_removed_file =>  $n_removed_file,);

isa_ok($n_remover, 'Genome::Model::Tools::Fasta::RemoveN');
