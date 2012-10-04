#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Temp;
require File::Compare;
use Test::More;

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx';
use_ok('Genome::Model::Tools::Sx::PhredReader') or die;
use_ok('Genome::Model::Tools::Sx::PhredWriter') or die;

is(Genome::Model::Tools::Sx::PhredReader->type, 'phred', 'type is phred');
is(Genome::Model::Tools::Sx::PhredWriter->type, 'phred', 'type is phred');

my $input_fasta = $dir.'/reader_writer.fasta';
ok(-s $input_fasta, 'Input fasta exists') or die;
my $input_qual = $input_fasta.'.qual';
ok(-s $input_qual, 'Input qual exists') or die;
my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
ok(-d $tmpdir, 'Created temp dir');
my $output_fasta = $tmpdir.'/fasta';
my $output_qual = $tmpdir.'/qual';

is_deeply([map { $_->property_name } Genome::Model::Tools::Sx::PhredReader->_file_properties], [qw/ file qual_file /], 'reader file properties');
is_deeply([map { $_->property_name } Genome::Model::Tools::Sx::PhredWriter->_file_properties], [qw/ file qual_file /], 'writer file properties');

my $failed_write = eval{ Genome::Model::Tools::Sx::PhredWriter->write(); };
diag($@);
ok(($@ && !$failed_write), 'Failed to write w/o fastqs');

# read/write w/ qual
my $reader = Genome::Model::Tools::Sx::PhredReader->create(
    file => $input_fasta,
    qual_file => $input_qual,
);
ok($reader, 'Create reader');

my $writer = Genome::Model::Tools::Sx::PhredWriter->create(
    file => $output_fasta,
    qual_file => $output_qual,
);
ok($writer, 'Create writer');

my $count = 0;
while ( my $seq = $reader->read ) {
    $count++;
    $writer->write($seq);
}
is($count, 10, 'Read/write 10 sequences');
ok($writer->flush, 'flush');
ok(-s $output_fasta, 'output fasta exists');
is(File::Compare::compare($input_fasta, $output_fasta), 0, 'In/out-put fastas match');
ok(-s $output_qual, 'output qual exists');
is(File::Compare::compare($input_qual, $output_qual), 0, 'In/out-put qauls match');

# read/write w/o qual
$reader = Genome::Model::Tools::Sx::PhredReader->create(
    file => $input_fasta,
);
ok($reader, 'Create reader w/o qual');

unlink $output_fasta;
$writer = Genome::Model::Tools::Sx::PhredWriter->create(
    file => $output_fasta,
);
ok($writer, 'Create writer w/o qual');

$count = 0;
while ( my $seq = $reader->read ) {
    $count++;
    $writer->write($seq);
}
is($count, 10, 'Read/write 10 sequences');
ok($writer->flush, 'flush');
ok(-s $output_fasta, 'output fasta exists');
is(File::Compare::compare($input_fasta, $output_fasta), 0, 'In/out-put fastas match');

# read fails
my $id_does_not_match = $dir.'/reader_writer.id_does_not_match.fasta.qual';
$reader = Genome::Model::Tools::Sx::PhredReader->create(
    file => $input_fasta, 
    qual_file => $id_does_not_match,
);
ok($reader, 'Create reader');
my $rv = eval{ $reader->read; };
diag($@);
ok((!$rv && $@ =~ /^Fasta and quality ids do not match:/), 'Failed when base and quals do not match');
my $quals_do_not_match = $dir.'/reader_writer.quals_do_not_match.fasta.qual';
$reader = Genome::Model::Tools::Sx::PhredReader->create(
    file => $input_fasta, 
    qual_file => $quals_do_not_match,
);
ok($reader, 'Create reader');
$rv = eval{ $reader->read; };
diag($@);
ok((!$rv && $@ =~ /^Number of qualities does not match/), 'Failed when id in fasta does not match qual');

#print $tmpdir; <STDIN>;
done_testing();
exit;

