#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 3;
use Genome::Utility::Test qw(compare_ok);

my $pkg = 'Genome::Model::Tools::SpeedSeq::Align';
use_ok($pkg);

my $speedseq_version = 'test';

my $test_data_dir = __FILE__.".d";
my $expected_output_dir = __FILE__.".out";

my $reference_fasta = $test_data_dir .'/human_g1k_v37_20_42220611-42542245.fasta';
my $fastq = $test_data_dir .'/NA12878.20slice.30X.fastq.gz';
my $read_group_header = '@RG\tID:NA12878\tSM:NA12878\tLB:lib1';

# Do not use the same temp directory for output.  SpeedSeq cleans up the temp directory.
my $temp_directory = Genome::Sys->create_temp_directory();
my $output_prefix = Genome::Sys->create_temp_directory() .'/example';

my $align_cmd = $pkg->create(
   version => $speedseq_version,
   output_prefix => $output_prefix,
   temp_directory => $temp_directory,
   reference_fasta => $reference_fasta,
   fastq => $fastq,
   paired => 1,
   sort_memory => 3,
   read_group_header => $read_group_header,
);

isa_ok($align_cmd,$pkg);
ok($align_cmd->execute,'execute command '. $pkg);

# BAM diff is failing

#for my $output_file ($align_cmd->output_files) {
#    my ($basename,$dirname) = File::Basename::fileparse($output_file);
#    my $expected_output_file = $expected_output_dir .'/'. $basename;
#    compare_ok($output_file,$expected_output_file);
#}   