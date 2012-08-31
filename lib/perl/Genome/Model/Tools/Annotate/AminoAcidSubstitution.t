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

#my $output = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Annotate-AminoAcidSubstitution/output";
my $output_fh = File::Temp->new(
    TEMPLATE => 'Genome-Model-Tools-Annotate-AminoAcidSubstitution-XXXXXX',
    DIR => "$ENV{GENOME_TEST_TEMP}/",
    CLEANUP => 1,
    UNLINK => 1,
);
my $output = $output_fh->filename;
$output_fh->close;

my $AminoAcidSubstitution = Genome::Model::Tools::Annotate::AminoAcidSubstitution->create(
    transcript => "ENST00000269305", 
    amino_acid_substitution => "S166C", 
    organism => "human", 
    output => $output
);

ok ($AminoAcidSubstitution, "successfully created command object");

my $gt_dir = "/gscmnt/200/medseq/biodb/shared/misc/annotation/" . $AminoAcidSubstitution->version;
ok(-d $gt_dir, "Gene/Transcript directory exists at $gt_dir");

my ($amino_acid_substitution) = $AminoAcidSubstitution->execute();
ok ($amino_acid_substitution, "got result from command execution");
ok (-e "$output.txt", "output file exists");

my @command = ["diff" , $expected , "$output.txt"];
my ($out) = &ipc_run(@command);
ok (! $out, "output matches expected");

sub ipc_run {
    my (@command) = @_;
    my ($in, $out, $err);
    IPC::Run::run(@command, \$in, \$out, \$err);
    
    return unless $out;
    return $out;	    
}
