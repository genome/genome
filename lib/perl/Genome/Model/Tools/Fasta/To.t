#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 8;
use File::Compare;
use File::Temp;

BEGIN {
        use_ok ('Genome::Model::Tools::Fasta::To::Fastq');
        use_ok ('Genome::Model::Tools::Fasta::To::Phd');
        use_ok ('Genome::Model::Tools::Fasta::To::Scf');
}

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fasta/To';
my $fasta_file = $dir .'/test.fasta';
my $expected_fastq = $dir .'/test.fastq.ori';

my $tmp_dir = File::Temp::tempdir(
    "FastaTo_XXXXXX", 
    TMPDIR => 1,
    CLEANUP => 1,
);

my %params = (
    fasta_file => $fasta_file,
    dir        => $tmp_dir,
);

my $to_fastq = Genome::Model::Tools::Fasta::To::Fastq->create(%params);

my $fastq_file = $to_fastq->out_file;

isa_ok($to_fastq,'Genome::Model::Tools::Fasta::To::Fastq');
ok($to_fastq->execute,'To::Fastq executes ok');
is(compare($fastq_file,$expected_fastq),0,'fastq files the same');
unlink $fastq_file || warn "Failed to remove fastq file $fastq_file";

my $to_phd = Genome::Model::Tools::Fasta::To::Phd->create(%params);
my @phd_list = $to_phd->list;
ok($to_phd->execute, 'To::Phd excutes ok');

my $to_scf = Genome::Model::Tools::Fasta::To::Scf->create(%params);
my @scf_list = $to_scf->list;
ok($to_scf->execute, 'To::Scf excutes ok');

my @list = map{$dir."/$_"}(@phd_list, @scf_list);
map{unlink}@list;
