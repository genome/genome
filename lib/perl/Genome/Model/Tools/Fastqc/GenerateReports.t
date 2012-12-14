#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 3;

use above 'Genome';

use_ok('Genome::Model::Tools::Fastqc::GenerateReports');
my $tmp_dir = Genome::Sys->create_temp_directory;
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fastqc-GenerateReports';

my @fastq_files = ($data_dir .'/s_2_1_sequence.txt',  $data_dir .'/s_2_2_sequence.txt');
my $fastq_files = join(',',@fastq_files);
my $fastqc = Genome::Model::Tools::Fastqc::GenerateReports->create(
    fastq_files      => $fastq_files,
    report_directory => $tmp_dir,
    use_version      => '0.4.3',
);
isa_ok($fastqc,'Genome::Model::Tools::Fastqc::GenerateReports');
ok($fastqc->execute,'execute command '. $fastqc->command_name);

#TODO: Add file comparsion or another test to verify output is complete and correct
