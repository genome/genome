#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 4;

use above 'Genome';

use_ok('Genome::InstrumentData::Command::EstimateLibraryComplexity');

my $test_output_dir = File::Temp::tempdir('Genome-InstrumentData-Command-EstimateLibraryComplexity-XXXXX', DIR => "$ENV{GENOME_TEST_TEMP}", CLEANUP => 1);

my $output_file = join('/', $test_output_dir, 'result.txt');

my $instrument_data_id = 2851949190;

my $command = Genome::InstrumentData::Command::EstimateLibraryComplexity->create(
    instrument_data_id => $instrument_data_id,
    output_file => $output_file,
    picard_version => '1.21',
);

ok($command, 'Created library complexity estimation command.');
ok($command->execute, 'Executed library complexity estimation command.');

#Beware--this string contains tabs (\t)!
my $expected_result = <<EORESULT
## net.sf.picard.metrics.StringHeader
# net.sf.picard.sam.EstimateLibraryComplexity OUTPUT=XXXXX MIN_IDENTICAL_BASES=5 MAX_DIFF_RATE=0.03 MIN_MEAN_QUALITY=20 OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 TMP_DIR=XXXXX VALIDATION_STRINGENCY=SILENT    READ_NAME_REGEX=[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).* VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000
## net.sf.picard.metrics.StringHeader
# Started on: XXXXX

## METRICS CLASS	net.sf.picard.sam.DuplicationMetrics
LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
H_HY-03023-lib1	0	1966	0	0	9	9	0.004578	

## HISTOGRAM	java.lang.Integer
duplication_group_count	H_HY-03023-lib1
1	1948
2	9

EORESULT
;

my $diff = Genome::Sys->diff_file_vs_text($output_file, $expected_result);

my @diff_lines = split("\n", $diff);
@diff_lines = grep($_ !~ m/XXXXX/, @diff_lines);

ok(!scalar(@diff_lines), 'Output matched expected result')
    or diag("  diff:\n" . $diff);
