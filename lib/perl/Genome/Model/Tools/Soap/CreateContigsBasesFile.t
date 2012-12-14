#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

use_ok ('Genome::Model::Tools::Soap::CreateContigsBasesFile') or die;

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Soap/CreateContigsBasesFile';

ok(-d $data_dir, "Data dir exists");

my $test_file = $data_dir.'/TEST.scafSeq';
ok(-s $test_file, "Test scaffold file exists");

my $temp_dir = Genome::Sys->create_temp_directory();
ok(-d $temp_dir, "Temp test dir created");

#link test.scafSeq file
symlink( $data_dir.'/TEST.scafSeq', $temp_dir.'/TEST.scafSeq' );
ok( -s $temp_dir.'/TEST.scafSeq', "Linked Test.scafSeq file");

my $create = Genome::Model::Tools::Soap::CreateContigsBasesFile->create(
    assembly_directory => $temp_dir,
    min_contig_length => 50,
    );
ok($create, "Created gmt soap create-contigs-fasta-file");
$create->dump_status_messages(1);
ok( ($create->execute) == 1, "Create executed successfully");

#check
my $file_name = 'contigs.bases';
ok(-s $data_dir."/$file_name", "Data dir contigs.fasta file exists");
ok(-s $temp_dir."/edit_dir/$file_name", "Created contigs fasta file");

#compare output file
ok(File::Compare::compare("$data_dir/$file_name", "$temp_dir/edit_dir/$file_name") == 0, "Output files match");

#<STDIN>;

done_testing();
