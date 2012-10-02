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

use_ok('Genome::Model::Tools::Sx::SffReader') or die;

is(Genome::Model::Tools::Sx::SffReader->type, 'sff', 'type is sff');

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
ok(-d $tmpdir, 'Created temp dir');
my $fastq = $tmpdir.'/out.fastq';

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/';
my $sff = $dir.'/reader.sff';
ok(-s $sff, 'sff exists') or die;
my $example_fastq = $dir.'/reader.sff.fastq';
ok(-s $example_fastq, 'example fastq exists') or die;

my $cmd = "gmt sx base -input $sff -output $fastq";
my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };

is(File::Compare::compare($example_fastq, $fastq), 0, 'fastq files match');

#print "gvimdiff $fastq $example_fastq\n"; <STDIN>;
done_testing();
exit;

