#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 5;

use_ok('Genome::Model::Tools::Graph::MutationCountBarPlot');

my $input_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Graph-MutationCountBarPlot";

# The input file
my $input_file_path = "$input_dir/input.txt";

# Temporary output files
my $output_data = Genome::Sys->create_temp_file_path;
my $output_R = Genome::Sys->create_temp_file_path;
my $output_plot = Genome::Sys->create_temp_file_path;

# Expected files
my $expected_data = "$input_dir/expected.dat";
my $expected_R = "$input_dir/expected.R";
my $expected_plot = "$input_dir/expected.png";

# Create a command
my $command= Genome::Model::Tools::Graph::MutationCountBarPlot->create(
    input_file_path => $input_file_path,
    output_file_path => $output_plot,
);

ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_plot, "output file created");

# The files will have a timestamp that will differ. Ignore this but check the rest.
my $diff = Genome::Sys->diff_file_vs_file($output_plot, $expected_plot);
ok($diff, 'output matched expected result')
