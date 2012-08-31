#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Temp;
require File::Compare;
use Test::More;

use_ok('Genome::Model::Tools::Sx::BedWriter') or die;

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
ok(-d $tmpdir, 'temp dir');
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/';
my $fasta = $dir.'/bed_writer.fasta';
ok(-s $fasta, 'fasta exists') or die;
my $example_bed_file = $dir.'/bed_writer.v2.bed';
ok(-s $example_bed_file, 'example bed file exists');

my $reader = Genome::Model::Tools::Sx::PhredReader->create(
    file => $fasta,
);
ok($reader, 'fasta reader');
my $bed_file = $tmpdir.'/bed';
my $writer = Genome::Model::Tools::Sx::BedWriter->create(
    file => $bed_file,
);
ok($writer, 'bed writer');
my $count = 0;
while ( my $fasta = $reader->read ) {
    $writer->write($fasta) or next;
    $count++;
}
is($count, 2, 'Read/write 2 beds');
is(File::Compare::compare($bed_file, $example_bed_file), 0, 'Bed file matches');

#print "$tmpdir\n"; <STDIN>;
done_testing();
exit;

