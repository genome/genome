#!/usr/bin/env perl

use above 'Genome';
use Test::More;
use File::Temp qw(tempdir);
use File::Slurp qw/write_file/;
use Data::Dumper;

use strict;
use warnings;

my $pkg = 'Genome::File::Vep::Command::ListGenes';
use_ok($pkg);

my $tmpdir = tempdir(CLEANUP => 1);
my $data = <<EOS
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	Extra
HELLO1_1_10_A_T	1:10	T	GENE1	TS1	Transcript	NON_SYNONYMOUS_CODING	3	3	1	D/V	gAt/gTt	-	EX1=e1;HGNC=HELLO1
HELLO1_1_11_A_T	1:10	T	GENE1	TS1	Transcript	NON_SYNONYMOUS_CODING	3	3	1	D/V	gAt/gTt	-	EX1=e1;HGNC=HELLO1
HELLO2_1_20_G_A	1:20	A	GENE2	TS2	Transcript	STOP_GAINED	6	6	2	W/*	tgG/tgA	-	HGNC=HELLO2
HELLO2_1_21_G_A	1:20	A	GENE2	TS2	Transcript	STOP_GAINED	6	6	2	W/*	tgG/tgA	-	HGNC=HELLO2
HELLO3_3_20_G_A	3:20	A	GENE3	TS3	Transcript	STOP_GAINED	6	6	2	W/*	tgG/tgA	-	HGNC=HELLO2
HELLO4_4_20_G_A	4:20	A	GENE4	TS4	Transcript	STOP_GAINED	6	6	2	W/*	tgG/tgA	-	HGNC=HELLO2;PolyPhen=possibly_damaging(0.6);Condel=neutral(0.1);SIFT=tolerated(0.2)
EOS
;

my $input_file_path = join("/", $tmpdir, "input.vep");
my $input_file = write_file($input_file_path, $data);

my $cmd = $pkg->create(input_file => $input_file_path);
ok($cmd, "Created command");
ok($cmd->execute, "Executed command");
my @genes = $cmd->genes;
@genes = sort(@genes);
my @expected = map {"GENE$_"} 1..4;
is_deeply(\@genes, \@expected, "Got expected genes") or diag(
    "Expected: " . Dumper(\@expected) . "\n" .
    "Actual: " . Dumper(\@genes));

done_testing();
