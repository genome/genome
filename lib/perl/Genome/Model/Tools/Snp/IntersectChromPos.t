#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 9;
use File::Compare;

use above 'Genome';

BEGIN {
    use_ok('Genome::Model::Tools::Snp::IntersectChromPos');
};
    
my $test_input_dir  = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Snp-IntersectChromPos/';

my $file1 = $test_input_dir . 'file1.in';
my $file2 = $test_input_dir . 'file2.in';

my $expected_intersect_output = $test_input_dir . 'intersect.expected';
my $expected_f1_only_output = $test_input_dir . 'f1_only.expected';
my $expected_f2_only_output = $test_input_dir . 'f2_only.expected';

my $test_output_dir = File::Temp::tempdir('Genome-Model-Tools-Snp-IntersectChromPos-XXXXX', CLEANUP => 1, TMPDIR => 1);
$test_output_dir .= '/';

my $intersect_output = $test_output_dir . 'intersect.out';
my $f1_only_output = $test_output_dir . 'f1_only.out';
my $f2_only_output = $test_output_dir . 'f2_only.out';

my $intersect_chrom_pos = Genome::Model::Tools::Snp::IntersectChromPos->create(
    file1 => $file1,
    file2 => $file2,
    intersect_output => $intersect_output,
    f1_only_output => $f1_only_output,
    f2_only_output => $f2_only_output
);

ok($intersect_chrom_pos, 'created IntersectChromPos object');
ok($intersect_chrom_pos->execute(), 'executed IntersectChromPos object');

ok(-s $intersect_output, 'intersection output generated');
is(compare($intersect_output, $expected_intersect_output), 0, 'intersection output matched expected result');

ok(-s $f1_only_output, 'file1 output generated');
is(compare($f1_only_output, $expected_f1_only_output), 0, 'file1 output matched expected result');

ok(-s $f2_only_output, 'file2 output generated');
is(compare($f2_only_output, $expected_f2_only_output), 0, 'file2 output matched expected result');
