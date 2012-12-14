#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Data::Dumper;
use Genome::Utility::Parser;

use Test::More tests => 5;
use Test::Differences;

# Use FindBin to provie a proper full path to the test files
# so we can run this test from anywhere
use FindBin qw($Bin);
my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Utility-Parser';
my $tsv_file = "$test_dir/test.tsv";
my $csv_file = "$test_dir/test.csv";

my @header = qw(build chromosome orientation start end sample allele1 allele2 comments);

# Tests using tab delimiter
my $tsv_parser = Genome::Utility::Parser->create(
                                                  file => $tsv_file,
                                                  separator => "\t",
                                             );
isa_ok($tsv_parser,'Genome::Utility::Parser');
is_deeply($tsv_parser->header_fields,\@header,'header parsed correctly for tsv');
my $line = 0;
my %tsv_data;
while (my $record = $tsv_parser->next) {
    $tsv_data{$line++} = $record;
}
$tsv_parser->close;

# Tests using comma delimiter
my $csv_parser = Genome::Utility::Parser->create(
                                                      file => $csv_file,
                                                  );
isa_ok($csv_parser,'Genome::Utility::Parser');
is_deeply($csv_parser->header_fields,\@header,'header parsed correctly for csv');

$line = 0;
my %csv_data;
while (my $record = $csv_parser->next) {
    $csv_data{$line++} = $record;
}
$csv_parser->close;

# Test equality of tab versus comma delimited
eq_or_diff(\%tsv_data,\%csv_data,'data produced by tab and comma delimited files');
