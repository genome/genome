#!/usr/bin/env genome-perl

use strict;
use warnings;

use File::Compare;
use Test::More tests => 4;

use Genome;

use_ok('Genome::Model::Tools::Fastq::ToFasta');

my $tmp_dir = Genome::Sys->create_temp_directory('Genome-Model-Tools-Fastq-ToFasta-'. Genome::Sys->username);
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fastq-To-Fasta';
my $fastq_file = $data_dir .'/s_6_1_sequence.txt';
my $expected_fasta_file = $data_dir .'/s_6_1_sequence.fa';

my $fastq_to_fasta = Genome::Model::Tools::Fastq::ToFasta->create(
   fastq_file => $fastq_file,
   fasta_file => $tmp_dir .'/s_6_1_sequence.fa',
);
isa_ok($fastq_to_fasta,'Genome::Model::Tools::Fastq::ToFasta');
ok($fastq_to_fasta->execute,'execute command '. $fastq_to_fasta->command_name);
ok(!compare($expected_fasta_file,$fastq_to_fasta->fasta_file),'expected fasta file '. $expected_fasta_file .' is identical to '. $fastq_to_fasta->fasta_file);
