#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More tests => 8;
use Text::Diff;
use File::Temp;


BEGIN {
    use_ok('Genome::Utility::SeqcleanReport::Reader');
    use_ok('Genome::Utility::SeqcleanReport::Writer');
}

my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Utility-SeqcleanReport';
my $file = "$test_dir/test.cln";

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $out_file = "$tmp_dir/out.cln";

my $reader = Genome::Utility::SeqcleanReport::Reader->create(
                                                             file => $file,
                                                         );
isa_ok($reader,'Genome::Utility::SeqcleanReport::Reader');
is($reader->separator,"\t",'separator');
is($reader->file,$file,'file accessor');

my $writer = Genome::Utility::SeqcleanReport::Writer->create(
                                                             file => $out_file,
                                                         );
isa_ok($writer,'Genome::Utility::SeqcleanReport::Writer');
is($writer->file,$out_file,'file accessor');
while (my $record = $reader->next) {
    $writer->write_record($record);
}
$writer->close;
$reader->close;

my $diff = `diff -b  $file $out_file`;
is($diff,'','Files are the same');
