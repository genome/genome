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

my $VERSION = 0; # Bump this each time tests data changes

my $test_dir = File::Spec->join(Genome::Utility::Test->data_dir($class), "v$VERSION");
diag "Test data located at $test_dir\n";

my $bam_file = File::Spec->join($test_dir, '134053361.bam');
ok(-s $bam_file, "Found input Bam file");

my $input_vcf = File::Spec->join($test_dir, 'testdata-indel.vcf');
ok(-s $input_vcf, "Found input Vcf file");

my $reference_sequence_build = Genome::Model::Build->get(106942997);
my $ref_fasta = $reference_sequence_build->full_consensus_path('fa');
diag "Reference sequence found at: $ref_fasta\n";
ok(-s $ref_fasta, "Found Reference Sequence fasta");


my $output_dir = Genome::Sys->create_temp_directory();

my $cmd = $class->create(
    make_realigned_bam => 1,
    input_vcf => $input_vcf,
    ref_fasta => $ref_fasta,
    bam_file => $bam_file,
    output_directory => $output_dir,
);
ok($cmd->execute(), "Successfully ran command");

my $expected_output_dir = File::Spec->join($test_dir, 'output-with-realigned-bam');
test_bam_files_are_identical();
test_txt_files_are_identical();

done_testing();

sub test_bam_files_are_identical {
    my $filename = 'all_events_realigned.merged.sorted.bam';
    my $expected = File::Spec->join($expected_output_dir, $filename);
    ok(-s $expected, "Found expected output bam file: $expected");
    my $found = File::Spec->join($output_dir, $filename);
    ok(-s $found, "Found resulting bam file: $found");
    my $diff = `diff <(samtools view $expected) <(samtools view $found)`;
    ok(!defined($diff), "No differences in Bam files.");
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
    push @files, glob(File::Spec->join($base, '*' ,'*.txt'));
    push @files, glob(File::Spec->join($base, '*', '*', '*.txt'));
    return @files;
}

