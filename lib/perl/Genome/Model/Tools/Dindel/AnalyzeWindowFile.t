#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);
use File::Spec;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

my $class = 'Genome::Model::Tools::Dindel::AnalyzeWindowFile';
use_ok($class) || die;

my $VERSION = 0; # Bump this each time tests data changes

my $test_dir = File::Spec->join(Genome::Utility::Test->data_dir($class), "v$VERSION");
diag "Test data located at $test_dir\n";

my $window_file = File::Spec->join($test_dir, 'inputs', 'dindel_window_file.txt');
ok(-s $window_file, "Found input window file") || die;

my $library_metrics_file = File::Spec->join($test_dir, 'inputs', 'cigar_generated_indels.libraries.txt');
ok(-s $library_metrics_file, "Found input library metrics file") || die;

my $input_bam = File::Spec->join($test_dir, 'inputs', '134053361-region-limited.bam');
ok(-s $input_bam, "Found input Bam file") || die;

my $reference_sequence_build = Genome::Model::Build->get(106942997);
my $ref_fasta = $reference_sequence_build->full_consensus_path('fa');
diag "Reference sequence found at: $ref_fasta\n";
ok(-s $ref_fasta, "Found Reference Sequence fasta") || die;

my $output_dir = Genome::Sys->create_temp_directory();

my $cmd = $class->create(
    window_file => $window_file,
    library_metrics_file => $library_metrics_file,
    input_bam => $input_bam,
    output_prefix => File::Spec->join($output_dir, 'results'),
    ref_fasta => $ref_fasta,
);
ok($cmd->execute(), "Successfully ran command");

my $expected_output_dir = File::Spec->join($test_dir, 'results');

my $expected_sam = File::Spec->join($test_dir, 'results', 'results_realigned.merged.sam');
test_bam_file_identical_to_sam_file($cmd->output_bam_file, $expected_sam);
done_testing();

sub test_bam_file_identical_to_sam_file {
    my ($bam_file, $sam_file) = @_;

    Genome::Model::Tools::Sam::BamToSam->execute(
        bam_file => $bam_file,
        include_headers => 1,
    );
    $bam_file =~ s/bam$/sam/;
    compare_ok($sam_file, $bam_file, "Output Bam file meets expectations.");
}
