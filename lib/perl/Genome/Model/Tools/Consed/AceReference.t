#!/usr/bin/env genome-perl

use strict;
use warnings;
use above 'Genome';
use Test::More tests => 6;

use_ok('Genome::Model::Tools::Consed::AceReference');
my $ace_file = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Consed-AceReference/AML_Validation_22_Trios-43114216_43116216-Ensembl.ace";
ok (-f $ace_file);

my ($info) = Genome::Model::Tools::Consed::AceReference->execute(ace_file => $ace_file, name_and_number => 1, no_stdout => 1);
ok ($info);
my $ace_reference = $info->result;
ok ($ace_reference);
my $Contig_number = $ace_reference->{Contig_number};
ok ($Contig_number eq "Contig1");
my $reseqid = $ace_reference->{reseqid};
ok ($reseqid eq "AML_Validation_22_Trios-43114216_43116216-Ensembl.c1");

