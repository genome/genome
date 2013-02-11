#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 7;
use File::Temp;
use File::Copy;

use above 'Genome';

use_ok('Genome::Model::Tools::Fasta::Split');
my $cleanup = 1;
my $test_fasta_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fasta/Chunk/test1.fa';
my $tmp_dir = File::Temp::tempdir("FastaSplit1_XXXXXX",CLEANUP=>1, TMPDIR => 1);
my $tmp_test_fasta_file = $tmp_dir.'/test1.fa';
copy $test_fasta_file, $tmp_test_fasta_file;

my $number_of_files = 3;
my $split_cmd = Genome::Model::Tools::Fasta::Split->create(
   number_of_files => $number_of_files,
   fasta_file => $tmp_test_fasta_file,
   output_directory => $tmp_dir,
);
isa_ok($split_cmd,'Genome::Model::Tools::Fasta::Split');
ok($split_cmd->execute,'execute command '. $split_cmd->command_name);
my $fasta_file_ref = $split_cmd->_split_fasta_files;
my @fasta_files = @{$fasta_file_ref};
is(scalar(@fasta_files),$number_of_files,'The number of files matches');

my $tmp_dir_2 = File::Temp::tempdir("FastaSplit2_XXXXXX",CLEANUP=>1, TMPDIR => 1);
my $tmp_test_fasta_file2 = $tmp_dir_2 . '/test1.fa';
copy $test_fasta_file, $tmp_test_fasta_file2;

my $split_cmd_2 = Genome::Model::Tools::Fasta::Split->create(
    min_sequence => 33360,
    fasta_file => $tmp_test_fasta_file2,
    output_directory => $tmp_dir_2,
);
isa_ok($split_cmd_2,'Genome::Model::Tools::Fasta::Split');
ok($split_cmd_2->execute,'execute command '. $split_cmd_2->command_name);
my $fasta_file_ref_2 = $split_cmd_2->_split_fasta_files;
my @fasta_files_2 = @{$fasta_file_ref_2};
is(scalar(@fasta_files_2),$number_of_files,'The number of files matches');
