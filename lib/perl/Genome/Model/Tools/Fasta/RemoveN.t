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

my $fasta_file = File::Spec->join(Genome::Config::get('test_inputs'), 'Genome-Model-Tools-Fasta', 'file.fasta');
my $n_removed_file = Genome::Sys->create_temp_file_path();

my $n_remover = Genome::Model::Tools::Fasta::RemoveN->create(fasta_file     =>  $fasta_file,
                                                             n_removed_file =>  $n_removed_file,);

isa_ok($n_remover, 'Genome::Model::Tools::Fasta::RemoveN');
