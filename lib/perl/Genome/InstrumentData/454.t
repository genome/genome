#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Compare;
use Test::More;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

use_ok('Genome::InstrumentData::454') or die;

my $sample = Genome::Sample->__define__(
    name => '__TEST_SAMPLE__',
);
ok($sample, 'define sample');
my $library = Genome::Library->__define__(
    name => '__TEST_LIBRARY__',
    sample => $sample,
);
ok($library, 'define library');

my $id = -1111;
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-InstrumentData-454';
my $id454 = Genome::InstrumentData::454->__define__(
    # This is info from a real 454 index region
    id => -1111,
    seq_id => -1111,
    region_id => -2222,
    region_number => 2,
    library => $library,
    run_name => 'R_2010_01_09_11_08_12_FLX08080418_Administrator_100737113',
    index_sequence => 'AACGGAGTC',
    read_count => 6437,
    sff_file => "$dir/2852582718.sff",
) or die "Unable to create mock 454 inst data";
my $sff_file = $id454->sff_file;
ok(-s $sff_file, 'sff file exists');
is($sff_file, $dir."/2852582718.sff", 'sff file named correctly');
is($id454->run_identifier, $id454->run_name.'.'.$id454->region_number.'-'.$id454->index_sequence, 'run identifier');
is($id454->read_count, 6437, 'read count');
is($id454->total_reads, 6437, 'total reads');

# test fasta, qual, fastq files and names w/ this real sff
my %types_methods = (
    fasta => 'dump_fasta_file',
    qual => 'dump_qual_file',
    fastq => 'dump_sanger_fastq_files',
);
for my $type ( keys %types_methods ) {
    my $method = $types_methods{$type};
    my ($file) = $id454->$method;
    ok($file, "dumped $type file");
    like($file, qr/$id(-output)?.$type$/, "$type file name matches: $file");
    my $example_file = "$dir/2852582718.$type";
    ok(-s $example_file, "example $type file exists: $example_file");
    is(File::Compare::compare($file, $example_file), 0, "$type file matches example file");
}

done_testing();
exit;

