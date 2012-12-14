#!/usr/bin/env genome-perl

use strict;
use warnings;

use File::Compare;
use Test::More tests => 4;

use above 'Genome';

use_ok('Genome::Model::Tools::Kmer::DuplicationRate');

my $tmp_dir = Genome::Sys->create_temp_directory('Genome-Model-Tools-Kmer-DuplicationRate-'. Genome::Sys->username);
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Kmer-DuplicationRate';
my $read_1_fastq = $data_dir .'/s_7_1_sequence.txt';
my $read_2_fastq = $data_dir .'/s_7_2_sequence.txt';
my $expected_file = $data_dir .'/s_7_occratio.txt';

my $dup_rate = Genome::Model::Tools::Kmer::DuplicationRate->create(
    fastq_files => $read_1_fastq .','. $read_2_fastq,
    output_file => $tmp_dir .'/s_7_occratio.txt',
);
isa_ok($dup_rate,'Genome::Model::Tools::Kmer::DuplicationRate');
ok($dup_rate->execute,'execute command '. $dup_rate->command_name);
ok(!compare($expected_file,$dup_rate->output_file),'expected output file '. $expected_file .' is identical to '. $dup_rate->output_file);
