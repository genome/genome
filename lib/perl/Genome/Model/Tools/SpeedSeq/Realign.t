#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 3;
use Genome::Utility::Test qw(compare_ok);

my $pkg = 'Genome::Model::Tools::SpeedSeq::Realign';
use_ok($pkg);

my $speedseq_version = 'test';

my $test_data_dir = __FILE__.".d";
my $expected_output_dir = __FILE__.".out";

my $reference_fasta = $test_data_dir .'/human_g1k_v37_20_42220611-42542245.fasta';
my $bam = $test_data_dir .'/NA12878.20slice.30X.bam';

# Do not use the same temp directory for output.  SpeedSeq cleans up the temp directory.
my $temp_directory = Genome::Sys->create_temp_directory();
my $output_prefix = Genome::Sys->create_temp_directory() .'/example';

my $realign_cmd = $pkg->create(
   version => $speedseq_version,
   output_prefix => $output_prefix,
   temp_directory => $temp_directory,
   reference_fasta => $reference_fasta,
   bams => $bam,
   sort_memory => 3,
);

isa_ok($realign_cmd,$pkg);
ok($realign_cmd->execute,'execute command '. $pkg);

# The BAM headers differ since realign uses a temp pipe as the input to the BWAMEM command
# Figure out a way to diff a BAM excluding the header, it must exist somewhere

#for my $output_file ($realign_cmd->output_files) {
#    my ($basename,$dirname) = File::Basename::fileparse($output_file);
#    my $expected_output_file = $expected_output_dir .'/'. $basename;
#    compare_ok($output_file,$expected_output_file);
#}   