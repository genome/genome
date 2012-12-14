#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Temp;
require File::Compare;
use Test::More;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

use_ok('Genome::Model::Tools::Sx::BamReader') or die;

is(Genome::Model::Tools::Sx::BamReader->type, 'bam', 'type is bam');

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
ok(-d $tmpdir, 'Created temp dir');
my $fasta = $tmpdir.'/out.fasta';
my $qual = $tmpdir.'/out.qual';

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/Sam/v1';
my $bam = $dir.'/rw.bam';
ok(-s $bam, 'bam exists') or die;
my $example_fasta = $dir.'/example.fasta';
ok(-s $example_fasta, 'example fasta exists') or die;
my $example_qual = $dir.'/example.qual';
ok(-s $example_qual, 'example qual exists') or die;

my $cmd = "gmt sx -input $bam -output file=$fasta:qual_file=$qual";
my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };

is(File::Compare::compare($example_fasta, $fasta), 0, 'fasta files match');
is(File::Compare::compare($example_qual, $qual), 0, 'qual files match');

#print "$tmpdir\n"; <STDIN>;
done_testing();
