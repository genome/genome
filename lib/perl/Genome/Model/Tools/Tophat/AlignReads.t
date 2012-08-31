#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;
use File::Compare;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}
else {
    plan tests => 8;
}

use_ok('Genome::Model::Tools::Tophat::AlignReads');

my $input_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Tophat-AlignReads';

my $expected_data_dir = $input_data_dir .'/expected_output';
my $expected_junctions_file = $expected_data_dir .'/junctions.bed';
my $expected_insertions_file = $expected_data_dir .'/insertions.bed';
my $expected_deletions_file = $expected_data_dir .'/deletions.bed';
my $expected_left_kept_reads_info = $expected_data_dir .'/left_kept_reads.info';
my $expected_right_kept_reads_info = $expected_data_dir .'/right_kept_reads.info';

my $tmp_dir = File::Temp::tempdir('Tophat-AlignReads-'.Genome::Sys->username.'-XXXX',DIR => "$ENV{GENOME_TEST_TEMP}",CLEANUP => 1);

my $read_1_fastq_file = $input_data_dir .'/random_1_1.fq';
my $read_2_fastq_file = $input_data_dir .'/random_1_2.fq';
my $reference_path = $input_data_dir .'/TTN.bowtie';

my $aligner = Genome::Model::Tools::Tophat::AlignReads->create(
   alignment_directory => $tmp_dir,
   read_1_fastq_list => $read_1_fastq_file,
   read_2_fastq_list => $read_2_fastq_file,
   reference_path => $reference_path,
   insert_size => '50',
   use_version => '1.3.0',
   bowtie_version => '0.12.7',
);

isa_ok($aligner,'Genome::Model::Tools::Tophat::AlignReads');
ok($aligner->execute,'execute command '. $aligner->command_name);
   
ok((compare($aligner->junctions_file,$expected_junctions_file) == 0),'junctions are identical');
ok((compare($aligner->insertions_file,$expected_insertions_file) == 0),'insertions are identical');
ok((compare($aligner->deletions_file,$expected_deletions_file) == 0),'deletions are identical');
ok((compare($aligner->left_kept_reads_info,$expected_left_kept_reads_info) == 0),'left_kept_reads.info is identical');
ok((compare($aligner->right_kept_reads_info,$expected_right_kept_reads_info) == 0),'right_kept_reads.info is identical');

exit;
