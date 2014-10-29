#!/usr/bin/env genome-perl

use strict;
use warnings;

use File::Temp;
use above "Genome";
use Genome::Model::Tools::Assembly::UnpackTraces;

use Test::More tests => 1;

my $temp_dir = File::Temp::tempdir (CLEANUP => 1);

ok(Genome::Model::Tools::Assembly::UnpackTraces->execute (
							  trace_file => $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Assembly-UnpackTraces/unpack_test.tgz',
							  clip_vector => 1,
							  clip_quality => 1,
							  data_out_dir => $temp_dir,
							  zip_files => 1,
							  )->result
);
