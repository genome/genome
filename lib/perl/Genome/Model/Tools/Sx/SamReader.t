#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Temp;
require File::Compare;
use Test::More;

use_ok('Genome::Model::Tools::Sx::SamReader') or die;

is(Genome::Model::Tools::Sx::SamReader->type, 'sam', 'type is sam');

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
ok(-d $tmpdir, 'Created temp dir');
my $fasta = $tmpdir.'/out.fasta';
my $qual = $tmpdir.'/out.qual';

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/Sam/v1';
my $sam = $dir.'/rw.sam';
ok(-s $sam, 'sam exists');
my $example_fasta = $dir.'/example.fasta';
ok(-s $example_fasta, 'example fasta exists');
my $example_qual = $dir.'/example.qual';
ok(-s $example_qual, 'example qual exists');

my $cmd = "gmt sx base -input $sam -output file=$fasta:qual_file=$qual";
my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };

is(File::Compare::compare($example_fasta, $fasta), 0, 'fasta files match');
is(File::Compare::compare($example_qual, $qual), 0, 'qual files match');

# fail
my $bad = $dir.'/bad1.sam';
my $reader = Genome::Model::Tools::Sx::SamReader->create(file => $bad);
ok($reader, 'create sam reader for bad file');
ok(!eval{$reader->read;}, 'read bad sam failed');
like($@, qr|^Length of sequence|, 'correct error for seq length diff');

#print "$tmpdir\n"; <STDIN>;
done_testing();
exit;

