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

my $class = 'Genome::Model::Tools::Dindel::GetCigarIndels';
use_ok($class);

my $VERSION = 1; # Bump this each time test data changes

my $test_dir = File::Spec->join(Genome::Utility::Test->data_dir($class), "v$VERSION");
diag "Test data located at $test_dir\n";

my $bam_file = File::Spec->join($test_dir, '134053361-region-limited.bam');
ok(-s $bam_file, "Found input Bam file");

my $reference_sequence_build = Genome::Model::Build->get(106942997);
my $ref_fasta = $reference_sequence_build->full_consensus_path('fa');
diag "Reference sequence found at: $ref_fasta\n";
ok(-s $ref_fasta, "Found Reference Sequence fasta");

my $output_dir = Genome::Sys->create_temp_directory();

my $cmd = $class->create(
    input_bam => $bam_file,
    ref_fasta => $ref_fasta,
    output_directory => $output_dir,
);
ok($cmd->execute(), "Successfully ran command");

my $expected_output_dir = File::Spec->join($test_dir, 'expected_output');
test_txt_files_are_identical();

done_testing();

sub test_txt_files_are_identical{
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
    return @files;
}
