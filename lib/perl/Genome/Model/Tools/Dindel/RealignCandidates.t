#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 5;
use Genome::Utility::Test;

my $class = 'Genome::Model::Tools::Dindel::RealignCandidates';
use_ok($class);

my $VERSION = 0; # Bump this each time test data changes

my $test_dir = File::Spec->join(Genome::Utility::Test->data_dir($class), "v$VERSION");
diag "Test data located at $test_dir\n";

my $input_dindel = File::Spec->join($test_dir, 'input.dindel');
ok(-s $input_dindel, 'Found input Dindel file');

my $reference_sequence_build = Genome::Model::Build->get(106942997);
my $ref_fasta = $reference_sequence_build->full_consensus_path('fa');
diag "Reference sequence found at: $ref_fasta\n";
ok(-s $ref_fasta, "Found Reference Sequence fasta");

my $output_prefix = Genome::Sys->create_temp_file_path();
my $output_file = "$output_prefix.variants.txt";

my $cmd = $class->create(
    ref_fasta => $ref_fasta,
    variant_file => $input_dindel,
    output_prefix => $output_prefix,
);  
ok($cmd->execute(), 'Successfully ran command');

my $expected_output_file = File::Spec->join($test_dir, 'output.variants.txt');
my $diff = `diff $output_file $expected_output_file`;
ok(!$diff, 'Output files are identical');
