#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use File::Temp;
use Genome::Utility::Test qw(compare_ok);
use Test::More tests => 6;

use Genome::Utility::PSL::Reader;
use Genome::Utility::PSL::Writer;

my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Utility-PSL';
my $file = "$test_dir/test.psl";

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $out_file = "$tmp_dir/out.psl";

my $reader = Genome::Utility::PSL::Reader->create(
                                                   file => $file,
                                               );
isa_ok($reader,'Genome::Utility::PSL::Reader');
is($reader->separator,"\t",'separator');
is($reader->file,$file,'file accessor');

my $writer = Genome::Utility::PSL::Writer->create(
                                               file => $out_file,
                                           );
isa_ok($writer,'Genome::Utility::PSL::Writer');
is($writer->file,$out_file,'file accessor');
while (my $record = $reader->next) {
    $writer->write_record($record);
}
$writer->close;
$reader->close;

compare_ok($file, $out_file, 'Files are the same');
