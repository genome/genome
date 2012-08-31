#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

my $archos = `uname -a`;
if ($archos !~ /64/) {
    plan skip_all => "Must run from 64-bit machine";
} else {
    if($ENV{UR_RUN_LONG_TESTS}) {
        plan tests => 6;
    } else {
        plan skip_all => 'This test takes about an hour to run and thus is skipped.  Use `ur test run --long` to enable.';
    }
}

my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-DetectVariants-Maq/';
my $test_working_dir = File::Temp::tempdir('DetectVariants-MaqXXXXX', DIR => "$ENV{GENOME_TEST_TEMP}/", CLEANUP => 1);

my $bam_input = $test_dir . '/alignments/whole_rmdup.map';

my $expected_dir = $test_dir . '/expected.v3/'; #Updated to .v3 after discussion about our BED standard [insertions start the base before; indels the first removed base]

my $ref_seq_build = Genome::Model::Build::ImportedReferenceSequence->get(type_name => 'imported reference sequence', name => 'NCBI-human-build36');
ok($ref_seq_build, 'Got a reference sequence build') or die('Test cannot continue without a reference sequence build');
is($ref_seq_build->name, 'NCBI-human-build36', 'Got expected reference for test case');

my $ref_seq_input = $ref_seq_build->full_consensus_path('bfa');
ok(Genome::Sys->check_for_path_existence($ref_seq_input), 'Got a reference binary FASTA') or die('Test cannot continue without a reference binary FASTA');

my $version = '0.7.1';
my $snv_parameters = my $indel_parameters = '-q 1';

my $command = Genome::Model::Tools::DetectVariants::Maq->create(
    reference_sequence_input => $ref_seq_input,
    aligned_reads_input => $bam_input,
    version => $version,
    snv_params => $snv_parameters,
    indel_params => $indel_parameters,
    detect_snvs => 1,
    detect_indels => 1,
    output_directory => $test_working_dir,
);
ok($command, 'Created `gmt detect-variants maq` command');
ok($command->execute, 'Executed `gmt detect-variants maq` command');

my $diff_cmd = sprintf('diff -r -q %s %s', $test_working_dir, $expected_dir);

my $diff = `$diff_cmd`;
is($diff, '', 'No differences in output from expected result from running maq for this version and parameters');
