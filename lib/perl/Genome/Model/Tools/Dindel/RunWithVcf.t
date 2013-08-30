#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);
use File::Spec;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{NO_LSF} = 1;
};

my $class = 'Genome::Model::Tools::Dindel::RunWithVcf';
use_ok($class);

my $VERSION = 1; # Bump this each time tests data changes

my $test_dir = File::Spec->join(Genome::Utility::Test->data_dir($class), "v$VERSION");
diag "Test data located at $test_dir\n";

my $bam_file = File::Spec->join($test_dir, 'inputs', '134053361-region-limited.bam');
ok(-s $bam_file, "Found input Bam file") || die;

my $input_vcf = File::Spec->join($test_dir, 'inputs', 'testdata-indel.vcf');
ok(-s $input_vcf, "Found input Vcf file")|| die;

my $reference_sequence_build = Genome::Model::Build->get(106942997);
my $ref_fasta = $reference_sequence_build->full_consensus_path('fa');
diag "Reference sequence found at: $ref_fasta\n";
ok(-s $ref_fasta, "Found Reference Sequence fasta") || die;


my $output_dir = Genome::Sys->create_temp_directory();

my $cmd = $class->create(
    input_vcf => $input_vcf,
    ref_fasta => $ref_fasta,
    bam_file => $bam_file,
    output_directory => $output_dir,
);
ok($cmd->execute(), "Successfully ran command");

my $found_bam_file = File::Spec->join($output_dir, 'results', 'result_1_realigned.merged.bam');

my $expected_output_dir = File::Spec->join($test_dir, 'results');
my $expected_bam_file = $found_bam_file;
$expected_bam_file =~ s/$output_dir/$expected_output_dir/;

test_bam_files_are_identical($expected_bam_file, $found_bam_file);
test_txt_files_are_identical();

done_testing();

sub bam_to_sam {
    my ($bam_file) = @_;

    my $sam_file = $bam_file;
    $sam_file =~ s/bam$/sam/;
    unlink($sam_file) if -s $sam_file;

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

sub test_txt_files_are_identical {
    my @expected_files = sort(get_text_files($expected_output_dir));

    for my $file (@expected_files) {
        my $found_file = $file;
        $found_file =~ s/$expected_output_dir/$output_dir/;
        compare_ok($file, $found_file, "$found_file identical to $file");
    }
}

sub get_text_files {
    my $base = shift;
    my @files = glob(File::Spec->join($base, '*.txt'));
    push @files, glob(File::Spec->join($base, '*', '*', '*.txt'));
    return @files;
}

