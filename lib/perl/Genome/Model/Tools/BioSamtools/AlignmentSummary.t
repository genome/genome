#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use File::Compare;

use above 'Genome';

if ($] < 5.010) {
  plan skip_all => "this test is only runnable on perl 5.10+"
}
plan tests => 11;

use File::Compare;
use above 'Genome';

use_ok('Genome::Model::Tools::BioSamtools');
use_ok('Genome::Model::Tools::BioSamtools::AlignmentSummary');


my $tmp_dir = File::Temp::tempdir('BioSamtools-AlignmentSummary-'.Genome::Sys->username.'-XXXX',CLEANUP => 1, TMPDIR => 1);

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BioSamtools/AlignmentSummary';

my $bam_file = $data_dir .'/test.bam';
my $regions_file = $data_dir .'/test_regions_zero_based_start.bed';
my $expected_output_file = $data_dir .'/alignment_summary_4.tsv';

my $as = Genome::Model::Tools::BioSamtools::AlignmentSummary->create(
    output_file => $tmp_dir .'/alignment_summary.tsv',
    bam_file => $bam_file,
    bed_file => $regions_file,
    wingspan => 0,
);
isa_ok($as,'Genome::Model::Tools::BioSamtools::AlignmentSummary');
ok($as->execute,'execute AlignmentSummary command '. $as->command_name);

ok(!compare($expected_output_file,$as->output_file),'expected output file '. $expected_output_file .' is identical to '. $as->output_file);

my $wingspan_bam_file = $data_dir .'/wingspan.bam';
my $wingspan_bed_file = $data_dir .'/wingspan.bed';

my $w_0_file = Genome::Sys->create_temp_file_path('alignment_summary_w_0.tsv');
my $expected_w_0_file = $data_dir .'/alignment_summary_w_0.tsv';
my $w_0_as = Genome::Model::Tools::BioSamtools::AlignmentSummary->create(
    output_file => $w_0_file,
    bam_file => $wingspan_bam_file,
    bed_file => $wingspan_bed_file,
    wingspan => 0,      
);  
isa_ok($w_0_as,'Genome::Model::Tools::BioSamtools::AlignmentSummary');
ok($w_0_as->execute,'execute AlignmentSummary command '. $w_0_as->command_name);
ok(!compare($expected_w_0_file,$w_0_as->output_file),'expected output file '. $expected_w_0_file .' is identical to '. $w_0_as->output_file);

my $w_500_file = Genome::Sys->create_temp_file_path('alignment_summary_w_500.tsv');
my $expected_w_500_file = $data_dir .'/alignment_summary_w_500.tsv';
my $w_500_as = Genome::Model::Tools::BioSamtools::AlignmentSummary->create(
    output_file => $w_500_file,
    bam_file => $wingspan_bam_file,
    bed_file => $wingspan_bed_file,
    wingspan => 500,
);
isa_ok($w_500_as,'Genome::Model::Tools::BioSamtools::AlignmentSummary');
ok($w_500_as->execute,'execute AlignmentSummary command '. $w_500_as->command_name);
ok(!compare($expected_w_500_file,$w_500_as->output_file),'expected output file '. $expected_w_500_file .' is identical to '. $w_500_as->output_file);
