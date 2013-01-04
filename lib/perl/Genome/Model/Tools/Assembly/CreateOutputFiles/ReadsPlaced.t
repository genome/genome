#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::Assembly::CreateOutputFiles::ReadsPlaced' ) or die;

my $data_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Assembly-CreateOutputFiles";
ok(-d $data_dir, "Found data directory: $data_dir") or die;

#make test dir
my $temp_dir = Genome::Sys->create_temp_directory();
Genome::Sys->create_directory( $temp_dir.'/edit_dir' );

#copy test input files to temp dir
for my $file ('gap.txt', 'contigs.bases', 'readinfo.txt') {
    ok( -e $data_dir.'/edit_dir/'.$file, "Test dir $file exists" );
    ok( File::Copy::copy( $data_dir.'/edit_dir/'.$file, $temp_dir.'/edit_dir/'.$file ) == 1, "Copied $file to temp_dir")
}

#run tool
my $ec = system("chdir $temp_dir; gmt assembly create-output-files reads-placed --directory $temp_dir");
ok($ec == 0, "Command ran successfully");

#checkout files
ok( -s $data_dir.'/edit_dir/reads.placed', "Data dir reads.placed file exists" );
ok( -s $temp_dir.'/edit_dir/reads.placed', "New reads.placed file created" );
ok( File::Compare::compare( $data_dir.'/edit_dir/reads.placed', $temp_dir.'/edit_dir/reads.placed' ) == 0, "Files match" );

done_testing();
