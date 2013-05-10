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
my $fastq = $tmpdir.'/out.fastq';

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/Sam/v2';
my $bam = $dir.'/test.bam';
ok(-s $bam, 'bam exists') or die;
my $example_fastq = $dir.'/example.fastq';
ok(-s $example_fastq, 'example fastq exists') or die;

my $cmd = "gmt sx -input $bam:cnt=2 -output file=$fastq";
my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };

is(File::Compare::compare($example_fastq, $fastq), 0, 'fastq files match');

#print "$tmpdir\n"; <STDIN>;
done_testing();
