#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 4;

use File::Temp qw();

use_ok('Genome::Model::Tools::CopyNumber::PlotSegments');

my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-CopyNumber-PlotSegments";
my $output_file = Genome::Sys->create_temp_file_path;
my $input_dir = "$test_dir";
my $input_file = "$input_dir/input1.seg.cbs";

my $output_fh = File::Temp->new();
my $output_pdf = $output_fh->filename();
$output_fh->close();
diag $output_pdf;

my $command= Genome::Model::Tools::CopyNumber::PlotSegments->create(
	segment_files => $input_file,
	genome_build => 37,
	output_pdf => $output_pdf,
);

ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_pdf, "output file created");
