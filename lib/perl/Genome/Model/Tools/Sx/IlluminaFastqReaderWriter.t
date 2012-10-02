#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Temp;
require File::Compare;
use Test::More;

#< Use >#
use_ok('Genome::Model::Tools::Sx::IlluminaFastqReader') or die;
use_ok('Genome::Model::Tools::Sx::IlluminaFastqWriter') or die;

#< Type >#
is(Genome::Model::Tools::Sx::IlluminaFastqReader->type, 'illumina', 'type is illumina');
is(Genome::Model::Tools::Sx::IlluminaFastqWriter->type, 'illumina', 'type is illumina');

my $failed_write = eval{ Genome::Model::Tools::Sx::IlluminaFastqWriter->write(); };
diag($@);
ok(($@ && !$failed_write), 'Failed to write w/o fastqs');

#< Files >#
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/';
my $example_fastq = $dir.'/reader_writer.fastq';
ok(-s $example_fastq, 'example fastq exists') or die;

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
ok(-d $tmpdir, 'Created temp dir');
my $out_fastq = $tmpdir.'/out.fastq';

#< Read illumina, write illumina >#
note('Read illumina, converts to sanger, writer illumina, converts back');
my $reader = Genome::Model::Tools::Sx::IlluminaFastqReader->create(
    file => $example_fastq,
);
ok($reader, 'Create reader');
my $writer = Genome::Model::Tools::Sx::IlluminaFastqWriter->create(
    file => $out_fastq,
);
ok($writer, 'Create writer');
my $count = 0;
while ( my $fastq = $reader->read ) {
    $writer->write($fastq) or next;
    $count++;
}
is($count, 25, 'Read/write 25 fastq sets');
ok($writer->flush, 'flush');
is(File::Compare::compare($out_fastq, $example_fastq), 0, 'files match');

#print "$tmpdir\n"; <STDIN>;
done_testing();
exit;

