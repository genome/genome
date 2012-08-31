#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

use_ok ('Genome::Model::Tools::Soap::CreateSupercontigsFastaFile');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Soap/CreateSupercontigsFastaFile';
ok(-d $data_dir, "Data dir exists") or die;
ok(-s $data_dir.'/TEST.scafSeq', "Data dir test file exists") or die;
ok(-s $data_dir.'/supercontigs.fasta', "Data dir supercontigs.fasta file exists") or die;

my $temp_dir = Genome::Sys->create_temp_directory();
ok(-d $temp_dir, "Temp test dir created") or die;

symlink( $data_dir.'/TEST.scafSeq', $temp_dir.'/TEST.scafSeq' );
ok( -s $temp_dir.'/TEST.scafSeq', "Linked TEST.scafSeq file");

#mkdir $temp_dir.'/edit_dir';
#ok(-d $temp_dir.'/edit_dir', "Created temp test dir edit_dir");

my $create = Genome::Model::Tools::Soap::CreateSupercontigsFastaFile->create(
    assembly_directory => $temp_dir,
    min_contig_length => 50,
    );
ok($create, "Created gmt soap create-supercontigs-fasta-file") or die;
ok(($create->execute) == 1, "Executed command") or die;
ok(-s $temp_dir.'/edit_dir/supercontigs.fasta', "Created supercontigs.fasta") or die;

#compare output files
ok(File::Compare::compare($data_dir.'/supercontigs.fasta', $temp_dir.'/edit_dir/supercontigs.fasta') ==0, "Output files match"); 

#<STDIN>;

done_testing();

exit;
