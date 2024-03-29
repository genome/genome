#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
use File::Compare;

$ENV{UR_DBI_NO_COMMIT} = 1;
$ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;

BEGIN {
    my $archos = `uname -a`;
    if ($archos !~ /64/) {
        plan skip_all => "Must run from 64-bit machine";
    } else {
        plan tests => 16;
    }
};

use_ok( 'Genome::Model::Tools::Somatic::HighConfidence');

my $test_input_dir  = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-Somatic-HighConfidence/';

my $sniper_file     = $test_input_dir . 'sniper.in';
my $tumor_bam_file  = $test_input_dir . 'tumor.tiny.bam';

my $map50_somatic50_expected = $test_input_dir . 'map50_somatic50.expected';
my $map10_somatic50_expected = $test_input_dir . 'map10_somatic50.expected';
my $map10_somatic80_expected = $test_input_dir . 'map10_somatic80.expected';

my $test_output_dir = File::Temp::tempdir('Genome-Model-Tools-Somatic-HighConfidence-XXXXX', CLEANUP => 1, TMPDIR => 1);
$test_output_dir .= '/';

my $output_file_default     = $test_output_dir . 'high_confidence-default.out';
my $output_file_m50_s50     = $test_output_dir . 'high_confidence-m50_s50.out';
my $output_file_m10_s50     = $test_output_dir . 'high_confidence-m10_s50.out';
my $output_file_m10_s80     = $test_output_dir . 'high_confidence-m10_s80.out';

my $high_confidence_default = Genome::Model::Tools::Somatic::HighConfidence->create(
    sniper_file         => $sniper_file,
    tumor_bam_file      => $tumor_bam_file,
    output_file         => $output_file_default,
);

ok($high_confidence_default, 'created HighConfidence object (default mapping & somatic quality)');
ok($high_confidence_default->execute(), 'executed HighConfidence object');

ok(-e $output_file_default, 'generated an (possibly empty) output file');


my $high_confidence_m50_s50 = Genome::Model::Tools::Somatic::HighConfidence->create(
    sniper_file         => $sniper_file,
    tumor_bam_file      => $tumor_bam_file,
    output_file         => $output_file_m50_s50,
    min_mapping_quality => 50,
    min_somatic_quality => 50,
);

ok($high_confidence_m50_s50, 'created HighConfidence object (quality minima: mapping 50, quality 50)');
ok($high_confidence_m50_s50->execute(), 'executed HighConfidence object');

ok(-s $output_file_m50_s50, 'generated an output file');
is(compare($output_file_m50_s50, $map50_somatic50_expected), 0, 'output matched expected output');

my $high_confidence_m10_s50 = Genome::Model::Tools::Somatic::HighConfidence->create(
    sniper_file         => $sniper_file,
    tumor_bam_file      => $tumor_bam_file,
    output_file         => $output_file_m10_s50,
    min_mapping_quality => 10,
    min_somatic_quality => 50,
);

ok($high_confidence_m10_s50, 'created HighConfidence object (quality minima: mapping 10, quality 50)');
ok($high_confidence_m10_s50->execute(), 'executed HighConfidence object');

ok(-s $output_file_m10_s50, 'generated an output file');
is(compare($output_file_m10_s50, $map10_somatic50_expected), 0, 'output matched expected output');

my $high_confidence_m10_s80 = Genome::Model::Tools::Somatic::HighConfidence->create(
    sniper_file         => $sniper_file,
    tumor_bam_file      => $tumor_bam_file,
    output_file         => $output_file_m10_s80,
    min_mapping_quality => 10,
    min_somatic_quality => 80,
);

ok($high_confidence_m10_s80, 'created HighConfidence object (quality minima: mapping 10, quality 80)');
ok($high_confidence_m10_s80->execute(), 'executed HighConfidence object');

ok(-s $output_file_m10_s80, 'generated an output file');
is(compare($output_file_m10_s80, $map10_somatic80_expected), 0, 'output matched expected output');
