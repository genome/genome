#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::Assembly::CreateOutputFiles' ) or die;
my $data_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Assembly-CreateOutputFiles";
ok(-d $data_dir, "Found data directory: $data_dir") or die;

#create test dir
my $temp_dir = Genome::Sys->create_temp_directory();
Genome::Sys->create_directory( $temp_dir.'/edit_dir' );

#copy files
for my $file ('contigs.bases', 'gap.txt') {
    ok (-s $data_dir.'/edit_dir/'.$file, "Test $file exists");
    ok(File::Copy::copy($data_dir."/edit_dir/$file", $temp_dir.'/edit_dir'),"Copied $file to temp dir");
}

#run
my $ec = system("chdir $temp_dir; gmt assembly create-output-files supercontigs-fasta --directory $temp_dir");
ok($ec == 0, "Command ran successfully");

#compare output files
ok( -s $data_dir.'/edit_dir/supercontigs.fasta', "Test supercontigs.fasta file exists" );
ok( -s $temp_dir.'/edit_dir/supercontigs.fasta', "New supercontigs.fasta file created" );
ok( File::Compare::compare($data_dir.'/edit_dir/supercontigs.fasta', $temp_dir.'/edit_dir/supercontigs.fasta') == 0, "Output files match" );

done_testing();
