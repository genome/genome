#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Temp;
require File::Compare;
use Test::More;

use_ok('Genome::Model::Tools::Sx::FastqReader') or die;
use_ok('Genome::Model::Tools::Sx::FastqWriter') or die;

is(Genome::Model::Tools::Sx::FastqReader->type, 'sanger', 'type is sanger');
is(Genome::Model::Tools::Sx::FastqWriter->type, 'sanger', 'type is sanger');

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
ok(-d $tmpdir, 'Created temp dir');
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/';

my $collated_fastq = $dir.'/reader_writer.collated.fastq';
ok(-s $collated_fastq, 'Collated fastq exists') or die;
my $forward_fastq = $dir.'/reader_writer.forward.fastq';
ok(-s $forward_fastq, 'Forward fastq exists') or die;
my $reverse_fastq = $dir.'/reader_writer.reverse.fastq';
ok(-s $reverse_fastq, 'Reverse fastq exists') or die;

my $failed_write = eval{ Genome::Model::Tools::Sx::FastqWriter->write(); };
diag($@);
ok(($@ && !$failed_write), 'Failed to write w/o fastqs');

my $reader = Genome::Model::Tools::Sx::FastqReader->create(
    file => $forward_fastq,
);
ok($reader, 'Create reader');
my $out_fastq = $tmpdir.'/test.fastq';
my $writer = Genome::Model::Tools::Sx::FastqWriter->create(
    file => $out_fastq,
);
ok($writer, 'Create writer');
my $count = 0;
while ( my $fastqs = $reader->read ) {
    $writer->write($fastqs) or next;
    $count++;
}
is($count, 12, 'Read/write 12 fastq sets');
ok($writer->flush, 'flush');
is(File::Compare::compare($forward_fastq, $out_fastq), 0, 'files match');

my $validate_good_reader = Genome::Model::Tools::Sx::FastqReader->create(file =>  $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/FastqReaderWriter/good.fastq');
ok($validate_good_reader, 'create reader to validate file');
ok($validate_good_reader->validate, 'validated good fastq');

my $validate_bad_reader = Genome::Model::Tools::Sx::FastqReader->create(file =>  $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/FastqReaderWriter/bad.fastq');
ok($validate_bad_reader, 'create reader to validate file');
ok(!$validate_bad_reader->validate, 'failed to validate bad fastq');

my $rv = eval{ $writer->write({ id => 'SEQ', seq => 'AATTGGCC', qual => 'abcdefg', }); };
ok(!$rv, 'Failed to write when seq and qual are not the same length');
ok($writer->write({ id => 'SEQ', seq => '', qual => '', }), 'write w/o seq and qual');

#print "$tmpdir\n"; <STDIN>;
done_testing();
