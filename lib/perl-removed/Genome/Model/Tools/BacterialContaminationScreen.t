#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More tests => 5; 

BEGIN {
    use_ok('Genome::Model::Tools::BacterialContaminationScreen');
}

my $path = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-BacterialContaminationScreen";
my $input_file = "$path/bcs_input.txt";  
my $output_file = Genome::Sys->create_temp_file_path("bcs_output.txt"); 
my $output_expected_file = "$path/bcs_output_expected.txt";

my $batch_screen = Genome::Model::Tools::BacterialContaminationScreen->create(
    type=>'read',
    input_file => $input_file,
    output_file => $output_file, 
);
isa_ok($batch_screen,'Genome::Model::Tools::BacterialContaminationScreen');

my $r = $batch_screen->execute;
ok($r, "executed successfully");

my $diff = `diff $output_file $output_expected_file`;
is($diff, '', "output from $output_file matches expected from $output_expected_file")
    or diag($diff);

ok( ! Genome::Model::Tools::BacterialContaminationScreen->create(), "Create w/o type - failed as expected");
