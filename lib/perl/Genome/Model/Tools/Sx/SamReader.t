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
my $fastq = $tmpdir.'/out.fastq';

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/Sam/v2';
my $sam = $dir.'/test.noheader.sam';
ok(-s $sam, 'sam exists');
my $example_fastq = $dir.'/example.fastq';
ok(-s $example_fastq, 'example fastq exists');

my $cmd = "gmt sx -input $sam:cnt=2 -output file=$fastq";
my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };

is(File::Compare::compare($example_fastq, $fastq), 0, 'fastq files match');

# fail
my $bad = $dir.'/bad.mismatch-seq-qual-length.sam';
my $reader = Genome::Model::Tools::Sx::SamReader->create(file => $bad);
ok($reader, 'create sam reader for bad file');
ok(!eval{$reader->read;}, 'read bad sam failed');
like($@, qr|^Length of sequence|, 'correct error for seq length diff');

#print "$tmpdir\n"; <STDIN>;
done_testing();
