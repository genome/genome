#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Temp;
require File::Compare;
use Test::More;

use_ok('Genome::Model::Tools::Sx::Reader') or die;
use_ok('Genome::Model::Tools::Sx::Writer') or die;

my $failed_create = Genome::Model::Tools::Sx::Reader->create();
ok(!$failed_create, 'Failed to create w/ reader w/o config');
$failed_create = Genome::Model::Tools::Sx::Writer->create();
ok(!$failed_create, 'Failed to create w/ writer w/o config ');

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
ok(-d $tmpdir, 'Created temp dir');
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/';

my $collated_fastq = $dir.'/reader_writer.collated.fastq';
ok(-s $collated_fastq, 'Collated fastq exists') or die;
my $forward_fastq = $dir.'/reader_writer.forward.fastq';
ok(-s $forward_fastq, 'Forward fastq exists') or die;
my $reverse_fastq = $dir.'/reader_writer.reverse.fastq';
ok(-s $reverse_fastq, 'Reverse fastq exists') or die;

# reader fails
my $fail = Genome::Model::Tools::Sx::Reader->create(
    config => [ 'name=no_file' ],
);
ok(!$fail, 'Failed to create reader w/o file');
$fail = Genome::Model::Tools::Sx::Reader->create(
    config => [ 'file=no.type' ],
);
ok(!$fail, 'Failed to create reader w/o type and file that type cannot be determined');

# writer fails
my $file = $tmpdir.'/file.fastq';
$fail = Genome::Model::Tools::Sx::Writer->create(
    config => [ 'name=no_file' ],
);
ok(!$fail, 'Failed to create writer w/o file');
$fail = Genome::Model::Tools::Sx::Writer->create(
    config => [ 'no.type' ],
);
ok(!$fail, 'Failed to create writer w/o type and file that type cannot be determined');
$fail = Genome::Model::Tools::Sx::Writer->create(
    config => [ $file.':name=pair', $file.':name=fwd', ],
);
ok(!$fail, 'Failed to create writer w/ pair and fwd');
$fail = Genome::Model::Tools::Sx::Writer->create(
    config => [ $file.':name=pair', $file.':name=rev', ],
);
ok(!$fail, 'Failed to create writer w/ pair and rev');
$fail = Genome::Model::Tools::Sx::Writer->create(
    config => [ $file.':name=pair', $file.':name=other', ],
);
ok(!$fail, 'Failed to create writer w/ pair and other named writer');
$fail = Genome::Model::Tools::Sx::Writer->create(
    config => [ $file.':name=fwd', ],
);
ok(!$fail, 'Failed to create writer w/ fwd and not rev');
$fail = Genome::Model::Tools::Sx::Writer->create(
    config => [ $file.':name=rev', ],
);
ok(!$fail, 'Failed to create writer w/ rev and not fwd');

# reader: types for files
ok(!eval{ Genome::Model::Tools::Sx::Reader->_type_for_file() }, 'type for undef file failed: '.$@);
ok(!Genome::Model::Tools::Sx::Reader->_type_for_file('file.unknown'), 'type for unknown file failed');
ok(Genome::Model::Tools::Sx::Reader->_type_for_file('-'), 'default type for STDOUT is sanger');
is(Genome::Model::Tools::Sx::Reader->_type_for_file('file.fastq'), 'sanger', 'type for fastq');
is(Genome::Model::Tools::Sx::Reader->_type_for_file('file.fq'), 'sanger', 'type for fq');
is(Genome::Model::Tools::Sx::Reader->_type_for_file('file.fasta'), 'phred', 'type for fasta');
is(Genome::Model::Tools::Sx::Reader->_type_for_file('file.fa'), 'phred', 'type for fa');
is(Genome::Model::Tools::Sx::Reader->_type_for_file('file.fna'), 'phred', 'type for fna');
is(Genome::Model::Tools::Sx::Reader->_type_for_file('file.efasta'), 'ephred', 'type for efasta');
is(Genome::Model::Tools::Sx::Reader->_type_for_file('file.sam'), 'sam', 'type for sam');
is(Genome::Model::Tools::Sx::Reader->_type_for_file('file.bam'), 'bam', 'type for bam');
is(Genome::Model::Tools::Sx::Reader->_type_for_file('file.sff'), 'sff', 'type for sff');

# writer: types for files
ok(!eval{ Genome::Model::Tools::Sx::Writer->_type_for_file() }, 'type for undef file failed: '.$@);
ok(!Genome::Model::Tools::Sx::Writer->_type_for_file('file.unknown'), 'type for unknown file failed');
ok(Genome::Model::Tools::Sx::Writer->_type_for_file('-'), 'default type for STDOUT is sanger');
is(Genome::Model::Tools::Sx::Writer->_type_for_file('file.fastq'), 'sanger', 'type for fastq');
is(Genome::Model::Tools::Sx::Writer->_type_for_file('file.fq'), 'sanger', 'type for fq');
is(Genome::Model::Tools::Sx::Writer->_type_for_file('file.fasta'), 'phred', 'type for fasta');
is(Genome::Model::Tools::Sx::Writer->_type_for_file('file.fa'), 'phred', 'type for fa');
is(Genome::Model::Tools::Sx::Writer->_type_for_file('file.fna'), 'phred', 'type for fna');
is(Genome::Model::Tools::Sx::Writer->_type_for_file('file.bed'), 'bed', 'type for bed');
is(Genome::Model::Tools::Sx::Writer->_type_for_file('file.fasta.gz'), 'phred', 'type for gzipped fasta');

# class for type
ok(!eval{ Genome::Model::Tools::Sx::Reader->_reader_class_for_type() }, 'class for undef file failed: '.$@);
ok(!Genome::Model::Tools::Sx::Reader->_reader_class_for_type('unknown'), 'class for unknown file failed');
is(Genome::Model::Tools::Sx::Reader->_reader_class_for_type('sanger'), 'Genome::Model::Tools::Sx::FastqReader', 'class for sanger');
is(Genome::Model::Tools::Sx::Reader->_reader_class_for_type('illumina'), 'Genome::Model::Tools::Sx::IlluminaFastqReader', 'class for illumina');
is(Genome::Model::Tools::Sx::Reader->_reader_class_for_type('phred'), 'Genome::Model::Tools::Sx::PhredReader', 'class for phred');
is(Genome::Model::Tools::Sx::Reader->_reader_class_for_type('ephred'), 'Genome::Model::Tools::Sx::PhredEnhancedSeqReader', 'class for ephred');

ok(!eval{ Genome::Model::Tools::Sx::Writer->_writer_class_for_type() }, 'class for undef file failed: '.$@);
ok(!Genome::Model::Tools::Sx::Writer->_writer_class_for_type('unknown'), 'class for unknown file failed');
is(Genome::Model::Tools::Sx::Writer->_writer_class_for_type('sanger'), 'Genome::Model::Tools::Sx::FastqWriter', 'class for sanger');
is(Genome::Model::Tools::Sx::Writer->_writer_class_for_type('illumina'), 'Genome::Model::Tools::Sx::IlluminaFastqWriter', 'class for illumina');
is(Genome::Model::Tools::Sx::Writer->_writer_class_for_type('phred'), 'Genome::Model::Tools::Sx::PhredWriter', 'class for phred');
is(Genome::Model::Tools::Sx::Writer->_writer_class_for_type('bed'), 'Genome::Model::Tools::Sx::BedWriter', 'class for bed');

# read paired, write paired
my $reader = Genome::Model::Tools::Sx::Reader->create(
    config => [ $forward_fastq, $reverse_fastq ],
);
ok($reader, 'Create reader');
my $out_collated_fastq = $tmpdir.'/test2.collated.fastq';
my $writer = Genome::Model::Tools::Sx::Writer->create(
    config => [ $out_collated_fastq ],
);
ok($writer, 'Create writer');
is($writer->{_strategy}, 'write_to_all');
my $count = 0;
while ( my $seqs = $reader->read ) {
    $writer->write($seqs) or next;
    $count++;
}
is($count, 12, 'Read/write 12 fastq sets');
is(File::Compare::compare($collated_fastq, $out_collated_fastq), 0, 'Reverse in/output files match');

# read pair, write separate
$reader->delete;
$reader = Genome::Model::Tools::Sx::Reader->create(
    config => [ $collated_fastq.':type=sanger:cnt=2' ],
);
ok($reader, 'Create reader');
my $out_forward_fastq = $tmpdir.'/test3.forward.fastq';
my $out_reverse_fastq = $tmpdir.'/test3.reverse.fastq';
$writer = Genome::Model::Tools::Sx::Writer->create(
    config => [ 'name=fwd:type=sanger:file='.$out_forward_fastq, 'name=rev:'.$out_reverse_fastq ],
);
ok($writer, 'Create writer');
is($writer->{_strategy}, 'write_forward_and_reverse_separately');
$count = 0;
while ( my $seqs = $reader->read ) {
    $writer->write($seqs) or next;
    $count++;
}
is($count, 12, 'Read/write 12 fastq sets');
is(File::Compare::compare($out_forward_fastq, $forward_fastq), 0, 'Foward in/output files match');
is(File::Compare::compare($out_reverse_fastq, $reverse_fastq), 0, 'Reverse in/output files match');

# read pair, write fwd, rev and sing #
$reader->delete;
$reader = Genome::Model::Tools::Sx::Reader->create(
    config => [ 'cnt=2:'.$collated_fastq ],
);
ok($reader, 'Create reader');
$out_forward_fastq = $tmpdir.'/test4.forward.fastq';
$out_reverse_fastq = $tmpdir.'/test4.reverse.fastq';
my $out_sing_fastq = $tmpdir.'/test4.sing.fastq';
$writer = Genome::Model::Tools::Sx::Writer->create(
    config => [ $out_forward_fastq.':name=fwd', $out_reverse_fastq.':name=rev', $out_sing_fastq.':name=sing' ],
);
ok($writer, 'Create writer');
is($writer->{_strategy}, 'write_forward_reverse_and_singletons_separately');
$count = 0;
while ( my $seqs = $reader->read ) {
    if ( $count++ == 8 ) {
        $seqs = [ $seqs->[1] ];
    }
    $writer->write($seqs) or die;
}
is($count, 12, 'Read/write 12 fastq sets');
my $example4_forward_fastq = $dir.'/reader_writer.example_4.forward.fastq';
ok(-s $example4_forward_fastq, 'Forward fastq example 4 exists');
is(File::Compare::compare($out_forward_fastq, $example4_forward_fastq), 0, 'foward output file matches');
my $example4_reverse_fastq = $dir.'/reader_writer.example_4.reverse.fastq';
ok(-s $example4_reverse_fastq, 'Reverse fastq example 4 exists');
is(File::Compare::compare($out_reverse_fastq, $example4_reverse_fastq), 0, 'reverse output file matches');
my $example4_sing_fastq = $dir.'/reader_writer.example_4.sing.fastq';
ok(-s $example4_sing_fastq, 'Reverse fastq example 4 exists') or die;
is(File::Compare::compare($out_sing_fastq, $example4_sing_fastq), 0, 'singleton output file matches');

# read pair, write fwd and rev sep
$reader->delete;
$reader = Genome::Model::Tools::Sx::Reader->create(
    config => [ 'cnt=2:'.$collated_fastq ],
);
ok($reader, 'Create reader');
$out_forward_fastq = $tmpdir.'/test5.forward.fastq';
$out_reverse_fastq = $tmpdir.'/test5.reverse.fastq';
$writer = Genome::Model::Tools::Sx::Writer->create(
    config => [ $out_forward_fastq.':name=fwd', $out_reverse_fastq.':name=rev' ],
);
ok($writer, 'Create writer');
is($writer->{_strategy}, 'write_forward_and_reverse_separately');
$count = 0;
while ( my $seqs = $reader->read ) {
    $writer->write($seqs) or die;
}
is(File::Compare::compare($out_forward_fastq, $forward_fastq), 0, 'foward file matches');
is(File::Compare::compare($out_reverse_fastq, $reverse_fastq), 0, 'reverse file matches');

# write singletons only
$reader->delete;
$reader = Genome::Model::Tools::Sx::Reader->create(
    config => [ 'cnt=2:'.$collated_fastq ],
);
ok($reader, 'Create reader');
$out_sing_fastq = $tmpdir.'/test6.sing.fastq';
$writer = Genome::Model::Tools::Sx::Writer->create(
    config => [ $out_sing_fastq.':name=sing' ],
);
ok($writer, 'Create writer');
is($writer->{_strategy}, 'write_singletons_only');
$count = 0;
while ( my $seqs = $reader->read ) {
    next if $count++ != 8;
    $seqs = [ $seqs->[1] ];
    $writer->write($seqs) or die;
}
is(File::Compare::compare($out_sing_fastq, $example4_sing_fastq), 0, 'singleton output file matches');

# write to named writers
$reader = Genome::Model::Tools::Sx::Reader->create(
    config => [ 'cnt=2:'.$collated_fastq ],
);
ok($reader, 'Create reader');
$out_forward_fastq = $tmpdir.'/test7.FWD.fastq';
$out_reverse_fastq = $tmpdir.'/test7.REVERSE.fastq';
$out_sing_fastq = $tmpdir.'/test7.SING.fastq';
$writer = Genome::Model::Tools::Sx::Writer->create(
    config => [ $out_forward_fastq.':name=FWD', $out_reverse_fastq.':name=REV', $out_sing_fastq.':name=SING' ],
);
ok($writer, 'Create writer');
is($writer->{_strategy}, 'write_to_named_writers');
$count = 0;
while ( my $seqs = $reader->read ) {
    if ( $count++ == 8 ) {
        $seqs->[1]->{writer_name} = 'SING';
        $seqs = [ $seqs->[1] ];
    }
    else {
        $seqs->[0]->{writer_name} = 'FWD';
        $seqs->[1]->{writer_name} = 'REV';
    }
    $writer->write($seqs) or die;
}
is(File::Compare::compare($out_forward_fastq, $example4_forward_fastq), 0, 'foward output files match');
is(File::Compare::compare($out_reverse_fastq, $example4_reverse_fastq), 0, 'reverse output files match');
is(File::Compare::compare($out_sing_fastq, $example4_sing_fastq), 0, 'singleton output files match');

# read fasta, write fasta and bed
my $example_fasta = $dir.'/bed_writer.fasta';
$reader->delete;
ok(-s $example_fasta, 'example fasta');
$reader = Genome::Model::Tools::Sx::Reader->create(
    config => [ 'file='.$example_fasta ],
);
ok($reader, 'Create reader');
my $out_fasta = $tmpdir.'/test8.fasta';
my $out_bed = $tmpdir.'/test8.bed';
$writer = Genome::Model::Tools::Sx::Writer->create(
    config => [ $out_fasta, 'type=bed:file='.$out_bed ],
);
ok($writer, 'Create writer');
is($writer->{_strategy}, 'write_to_all');
$count = 0;
while ( my $seqs = $reader->read ) {
    $writer->write($seqs) or die;
}
is(File::Compare::compare($out_fasta, $example_fasta), 0, 'fasta output file matches');
my $example_bed = $dir.'/bed_writer.v2.bed';
ok(-s $example_bed, 'example bed');
is(File::Compare::compare($out_bed, $example_bed), 0, 'bed output file matches');

#print "$tmpdir\n"; <STDIN>;
done_testing();
