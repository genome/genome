#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test 'compare_ok';

my $class = 'Genome::Model::Tools::Dindel::RealignCandidates';
use_ok($class);

my $VERSION = 0; # Bump this each time test data changes

my $test_dir = File::Spec->join(Genome::Utility::Test->data_dir($class), "v$VERSION");
diag "Test data located at $test_dir\n";

my $input_dindel_file = File::Spec->join($test_dir, 'input_dindel_file.variants.txt');
ok(-s $input_dindel_file, 'Found input Dindel file') || die;

my $reference_sequence_build = Genome::Model::Build->get(106942997);
my $ref_fasta = $reference_sequence_build->full_consensus_path('fa');
diag "Reference sequence found at: $ref_fasta\n";
ok(-s $ref_fasta, "Found Reference Sequence fasta") || die;

my $output_directory = Genome::Sys->create_temp_directory();

my $cmd = $class->create(
    input_dindel_file => $input_dindel_file,
    ref_fasta => $ref_fasta,
    output_directory => $output_directory,
);
ok($cmd->execute(), 'Successfully ran command');

my $expected_output_file = $cmd->output_dindel_file;
$expected_output_file =~ s/$output_directory/$test_dir/;
compare_ok($expected_output_file, $cmd->output_dindel_file, 'Output matches expected output');

done_testing();
