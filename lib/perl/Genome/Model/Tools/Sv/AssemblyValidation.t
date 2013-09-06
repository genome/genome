#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
use Genome::Utility::Test qw(compare_ok);
use File::Temp qw(tempfile);
use File::Path qw(rmtree);

BEGIN {
    my $archos = `uname -a`;
    if ($archos !~ /64/) {
        plan skip_all => "Must run from 64-bit machine";
    } 
    else {
        plan tests => 17;
    }
};

use_ok( 'Genome::Model::Tools::Sv::AssemblyValidation');

my $test_input_dir  = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sv-AssemblyValidation/v2/';
my $normal_bam  = $test_input_dir . 'normal.bam';
my $sv_file     = $test_input_dir . 'sv.file';

my @file_names = qw(normal.out normal.cm_aln.out normal.bp_seq.out);
my @expected_files = map{$test_input_dir . $_}@file_names;

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-Sv-AssemblyValidation-XXXXX', 
    TMPDIR => 1,
    CLEANUP => 1,
);

my $ref_build = Genome::Model::Build::ImportedReferenceSequence->get(name => "NCBI-human-build36");
my $ref_fa = $ref_build->full_consensus_path('fa');

my @test_out_files = map{$tmp_dir.'/test.'.$_}@file_names;

my $sv_valid = Genome::Model::Tools::Sv::AssemblyValidation->create(
    sv_file   => $sv_file,
    bam_files => $normal_bam,
    output_file => $test_out_files[0],
    cm_aln_file => $test_out_files[1],
    breakpoint_seq_file => $test_out_files[2],
    reference_file => $ref_fa,
);

ok($sv_valid, 'created AssemblyValidation object');
ok($sv_valid->execute(), 'executed AssemblyValidation object OK');

for my $i (0..2) {
    ok(-s $test_out_files[$i], 'generated output file: '.$file_names[$i].' ok');
    compare_ok($expected_files[$i], $test_out_files[$i], 'output matched expected results: ' . $file_names[$i],
        filters => qr(^#Bams: .*));
}

$test_input_dir = $test_input_dir."chromosomeBeginTest.v1/";
$normal_bam = $test_input_dir.'normal.bam';
my $tumor_bam = $test_input_dir.'tumor.bam';
$sv_file = $test_input_dir.'input_calls';

@file_names = qw(output.csv output.cm output.fasta);
@expected_files = map{$test_input_dir . $_}@file_names;

@test_out_files = map{$tmp_dir.'/test.'.$_}@file_names;

$sv_valid = Genome::Model::Tools::Sv::AssemblyValidation->create(
    sv_file => $sv_file,
    bam_files => "$normal_bam,$tumor_bam",
    output_file => $test_out_files[0],
    cm_aln_file => $test_out_files[1],
    breakpoint_seq_file => $test_out_files[2],
    reference_file => $ref_fa,
);

ok($sv_valid, 'created AssemblyValidation object');
ok($sv_valid->execute(), 'executed AssemblyValidation object');

for my $i (0..2) {
    ok(-s $test_out_files[$i], 'generated output file: ' . $file_names[$i]);
    compare_ok($expected_files[$i], $test_out_files[$i], 'output matched expected results: ' . $file_names[$i],
        filters => qr(^#Bams: .*));
}
