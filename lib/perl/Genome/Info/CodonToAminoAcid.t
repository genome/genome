#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Info::CodonToAminoAcid');

my %single_letter = Genome::Info::CodonToAminoAcid->single_letter;
my %three_letter  = Genome::Info::CodonToAminoAcid->three_letter;

#test single letter
is($single_letter{TTG}, 'L',  'Amino acid single letter of condon TTG is L');
is($single_letter{ACA}, 'T',  'Amino acid single letter of condon ACA is T');
is($single_letter{TAG}, 'X',  'Amino acid single letter of condon TAG is X (stop codon)');
is($single_letter{'-AA'}, 'indel', 'Amino acid single letter of condon -AA is indel');
is($single_letter{'T+A'}, 'refseq allele',  'Amino acid single letter of condon T+A is refseq allele');
is(scalar keys %single_letter, 208, 'There should be 208 codons');

#test three letter
is($three_letter{TTG}, 'Leu',  'Amino acid three letter of condon TTG is Leu');
is($three_letter{ACA}, 'Thr',  'Amino acid three letter of condon ACA is Thr');
is($three_letter{TAG}, 'AMB',  'Amino acid three letter of condon TAG is AMB (stop codon)');
is($three_letter{'-AA'}, 'indel', 'Amino acid three letter of condon -AA is indel');
is($three_letter{'T+A'}, 'refseq allele',  'Amino acid three letter of condon T+A is refseq allele');
is(scalar keys %three_letter, 208, 'There should be 208 codons');

done_testing();
