#!/usr/bin/env genome-perl

use strict;
use warnings;
use above "Genome";
use IPC::Run;
use Test::More tests => 8;

use_ok('Genome::Model::Tools::Annotate::AminoAcidSubstitution');
my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Annotate-AminoAcidSubstitution";
ok (-d $test_dir, "test dir exists at $test_dir");
my $expected = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Annotate-AminoAcidSubstitution/expected.txt";
ok (-e $expected, "expected output file exists at $expected");

my $tmpdir = File::Temp::tempdir(
    TEMPLATE => 'Genome-Model-Tools-Annotate-AminoAcidSubstitution-XXXXXX',
    TEMPDIR => 1,
    CLEANUP => 1,
);
my $output = join('/', $tmpdir, 'output.txt');
my $output_txt = "$output.txt";
diag $output;

my $AminoAcidSubstitution = Genome::Model::Tools::Annotate::AminoAcidSubstitution->create(
    transcript => "ENST00000269305",
    amino_acid_substitution => "S166C",
    organism => "human",
    output => $output,
);

ok ($AminoAcidSubstitution, "successfully created command object");

my $gt_dir = "/gscmnt/200/medseq/biodb/shared/misc/annotation/" . $AminoAcidSubstitution->version;
ok(-d $gt_dir, "Gene/Transcript directory exists at $gt_dir");

my ($amino_acid_substitution) = $AminoAcidSubstitution->execute();
ok ($amino_acid_substitution, "got result from command execution");
ok (-e $output_txt, "output file exists");

my $diff = Genome::Sys->diff_file_vs_file($expected, $output_txt);
ok(!$diff, 'output matches expected') || diag $diff;
