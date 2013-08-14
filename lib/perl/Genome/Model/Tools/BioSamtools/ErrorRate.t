#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use File::Compare;
use Test::More;

use_ok('Genome::Model::Tools::BioSamtools::ErrorRate');

my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BioSamtools/ErrorRate/';
my $bam = $test_dir . 'test.bam';

my $expected_out = $test_dir.'v1/expected.out';
my $output_file  = Genome::Sys->create_temp_file_path('test.error_rate.tsv');

my $cmd =Genome::Model::Tools::BioSamtools::ErrorRate->create(
    bam_file    => $bam,
    output_file => $output_file,
    version     => 0.7,
);

ok($cmd, 'gmt bio-samtools error-rate command created');
ok($cmd->execute, 'gmt bio-samtools error-rate command completed successfully');

is(compare($output_file, $expected_out), 0, 'Output file '.$output_file. ' is same as expected file '.$expected_out);

done_testing();
