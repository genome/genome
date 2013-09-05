#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Model::Tools::Dindel::TestHelpers qw(
    get_test_dir
    get_ref_fasta
    compare_output_to_test_data
);

my $class = 'Genome::Model::Tools::Dindel::RealignCandidates';
use_ok($class);

my $VERSION = 0; # Bump this each time test data changes

my $ref_fasta = get_ref_fasta();

my $test_dir = get_test_dir($class, $VERSION);
my $input_dindel_file = File::Spec->join($test_dir, 'input_dindel_file.variants.txt');
ok(-s $input_dindel_file, 'Found input Dindel file') || die;

my $output_directory = Genome::Sys->create_temp_directory();

my $cmd = $class->create(
    input_dindel_file => $input_dindel_file,
    ref_fasta => $ref_fasta,
    output_directory => $output_directory,
);
ok($cmd->execute(), 'Successfully ran command');

compare_output_to_test_data($cmd->output_dindel_file, $output_directory, $test_dir);

done_testing();
