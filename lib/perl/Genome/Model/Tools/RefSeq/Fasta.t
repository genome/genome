#!/usr/bin/env genome-perl

use strict;
use warnings;
use above 'Genome';
use Test::More tests => 11;

use_ok('Genome::Model::Tools::RefSeq::Fasta');
my $refseq_fasta = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-RefSeq-Fasta/AML_Validation_22_Trios-43114216_43116216-Ensembl.c1.refseq.fasta";
ok (-f $refseq_fasta);

my ($info) = Genome::Model::Tools::RefSeq::Fasta->execute(refseq_fasta => $refseq_fasta, no_stdout => 1);
ok ($info);
my $ref_head = $info->result;
ok ($ref_head);

my $orientation = $ref_head->{orientation};
ok ($orientation eq "plus");

my $genomic_coord = $ref_head->{genomic_coord};
ok ($genomic_coord eq "43114215");

my $chromosome = $ref_head->{chromosome};
ok ($chromosome eq "20");

my $length = $ref_head->{length};
ok ($length eq "2001");

my $start = $ref_head->{start};
ok ($start eq "43114216");

my $stop = $ref_head->{stop};
ok ($stop eq "43116216");

my $name = $ref_head->{name};
ok ($name eq "AML_Validation_22_Trios-43114216_43116216-Ensembl.c1.refseq.fasta");
