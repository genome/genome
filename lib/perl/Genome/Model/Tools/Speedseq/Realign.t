#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 9;
use Genome::Utility::Test qw(compare_ok);

my $pkg = 'Genome::Model::Tools::Speedseq::Realign';
use_ok($pkg);

my $speedseq_version = 'test';

my $test_data_dir = __FILE__.".d";
my $expected_output_dir = __FILE__.".out";

my $reference_fasta = $test_data_dir .'/human_g1k_v37_20_42220611-42542245.fasta';
my $bam = $test_data_dir .'/NA12878.20slice.30X.bam';

# Do not use the same temp directory for output.  speedseq cleans up the temp directory.
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

# BAM diff
for my $output_file ($realign_cmd->output_files) {
    my ($basename,$dirname,$suffix) = File::Basename::fileparse($output_file,qw/\.bam \.bai/);  
    my $expected_output_file = $expected_output_dir .'/'. $basename . $suffix;
    if ($suffix eq '.bam') {
          my $cmp = Genome::Model::Tools::Sam::Compare->execute(
          file1 => $output_file,
          file2 => $expected_output_file,
          );
          ok ($cmp->result, 'compare BAMs'); 
    } elsif ($suffix eq '.bai') {
      ok(-s $output_file, 'BAM index file has size.');
    }                                              
}   