#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test 'compare_ok';
use Genome::Model::Tools::Dindel::TestHelpers qw(
    get_test_dir
    get_ref_fasta
    compare_output_to_test_data
);

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

my $class = 'Genome::Model::Tools::Dindel::AnalyzeWindowFile';
use_ok($class) || die;

my $VERSION = 0; # Bump this each time tests data changes

my $ref_fasta = get_ref_fasta();

my $test_dir = get_test_dir($class, $VERSION);
my $window_file = File::Spec->join($test_dir, 'inputs', 'dindel_window_file.txt');
ok(-s $window_file, "Found input window file") || die;

my $library_metrics_file = File::Spec->join($test_dir, 'inputs', 'cigar_generated_indels.libraries.txt');
ok(-s $library_metrics_file, "Found input library metrics file") || die;

my $input_bam = File::Spec->join($test_dir, 'inputs', '134053361-region-limited.bam');
ok(-s $input_bam, "Found input Bam file") || die;

my $output_directory = Genome::Sys->create_temp_directory();

my $cmd = $class->create(
    window_file => $window_file,
    library_metrics_file => $library_metrics_file,
    input_bam => $input_bam,
    ref_fasta => $ref_fasta,
    output_directory => $output_directory,
);
ok($cmd->execute(), "Successfully ran command");

(my $expected_bam = $cmd->output_bam) =~ s/$output_directory/$test_dir/;
test_bam_files_are_identical($expected_bam, $cmd->output_bam);
compare_output_to_test_data($cmd->output_glf, $output_directory, $test_dir);
compare_output_to_test_data($cmd->output_log, $output_directory, $test_dir);

done_testing();

sub bam_to_sam {
    my ($bam_file) = @_;

    my $sam_file = $bam_file;
    $sam_file =~ s/bam$/sam/;
    return $sam_file if -f $sam_file;

    Genome::Model::Tools::Sam::BamToSam->execute(
        bam_file => $bam_file,
        include_headers => 1,
        sam_file => $sam_file,
    );
    return $sam_file;
}

sub test_bam_files_are_identical {
    my ($bam1, $bam2) = @_;
    ok(-s $bam1, "Bam file exists: $bam1");
    ok(-s $bam2, "Bam file exists: $bam2");

    compare_ok(bam_to_sam($bam1), bam_to_sam($bam2), "No differences in Bam files.");
}
