#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 11;
use Genome::File::Vcf::Differ;
use Genome::Utility::Test qw(compare_ok);
use Genome::Test::Data qw(get_test_file);

my $pkg = 'Genome::Model::Tools::Speedseq::Sv';
use_ok($pkg);

my $speedseq_version = 'test';

my $expected_output_dir = __FILE__.".out";

my $reference_fasta = get_test_file('NA12878', 'human_g1k_v37_20_42220611-42542245.fasta');
my $bam = get_test_file('NA12878', 'NA12878.20slice.30X.aligned.bam');
my $split_bam = get_test_file('NA12878','NA12878.20slice.30X.splitters.bam');
my $discordant_bam = get_test_file('NA12878','NA12878.20slice.30X.discordants.bam');

my $temp_directory = Genome::Sys->create_temp_directory();
my $output_prefix = Genome::Sys->create_temp_directory() .'/example';
# Do not use the same temp directory for output.  speedseq cleans up the temp directory.

my $sv_cmd = $pkg->create(
   version => $speedseq_version,

   reference_fasta => $reference_fasta,
   full_bam_file => $bam,
   output_prefix => $output_prefix,
   split_read_bam_file => $split_bam,
   discordant_read_bam_file => $discordant_bam,
);
isa_ok($sv_cmd,$pkg);
ok($sv_cmd->execute,'execute command '. $pkg);

#VCF Difference
my $output_file = "$output_prefix.sv.vcf.gz";

my $expected_output_file = ("$expected_output_dir/Sv.t.out.sv.vcf.gz");

compare_ok($output_file, $expected_output_file);
compare_ok("$output_file.tbi", "$expected_output_file.tbi");

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~Third Test CNVnator Check~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my $output_prefix3 = Genome::Sys->create_temp_directory() .'/example.CNV';
# Do not use the same temp directory for output.  speedseq cleans up the temp directory.
my $expected_output_prefix3 = ("$expected_output_dir/Sv.t.out3");


my $sv_cmd3 = $pkg->create(
   version => $speedseq_version,
   temp_directory => $temp_directory,
   reference_fasta => $reference_fasta,
   full_bam_file => $bam,
   output_prefix => $output_prefix3,
   CNVnator_read_depth => 'true',
   split_read_bam_file => $split_bam,
   genotype_svtyper => 1,
   discordant_read_bam_file => $discordant_bam,
);
isa_ok($sv_cmd3,$pkg);
ok($sv_cmd3->execute,'execute command '. $pkg);

#VCF Difference
my $output_file3 = "$output_prefix3.sv.vcf.gz";

my $expected_output_file3 = ("$expected_output_dir/Sv.t.out3.sv.vcf.gz");

my $differ = Genome::File::Vcf::Differ->new($output_file3, $expected_output_file3);
    my $diff = $differ->diff;
    is($diff, undef, "Found No differences between $output_file3 and (expected) $expected_output_file3") ||
       diag $diff->to_string;

#compare_ok($output_file3, $expected_output_file3);
compare_ok("$output_file3.tbi", "$expected_output_file3.tbi");
compare_ok("$output_prefix3.sv.NA12878.20slice.30X.aligned.bam.readdepth.bed", "$expected_output_prefix3.sv.NA12878.20slice.30X.aligned.bam.readdepth.bed");
compare_ok("$output_prefix3.sv.NA12878.20slice.30X.aligned.bam.readdepth.txt", "$expected_output_prefix3.sv.NA12878.20slice.30X.aligned.bam.readdepth.txt");

