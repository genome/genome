#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::Assembly::CreateOutputFiles::FastaAndQualFromVelvetAfg' ) or die;

my $data_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Assembly-CreateOutputFiles2";
ok(-d $data_dir, "Found data directory: $data_dir") or die;

#test afg file
my $afg_file = $data_dir.'/velvet_asm.afg';
ok(-s $afg_file, "Test afg file exists") or die;

#make temp dir
my $temp_dir = Genome::Sys->create_temp_directory();
Genome::Sys->create_directory( $temp_dir.'/edit_dir' );

#link afg file in tmp dir
symlink($data_dir.'/velvet_asm.afg', $temp_dir.'/velvet_asm.afg') or die; 
ok (-s $temp_dir.'/velvet_asm.afg', "Linked afg file in tmp dir") or die;

my $ec = system("chdir $temp_dir; gmt assembly create-output-files fasta-and-qual-from-velvet-afg --directory $temp_dir");
ok($ec == 0, "Command ran successfully") or die;

for my $file ('contigs.bases', 'contigs.quals') {
    ok( -s $data_dir."/edit_dir/$file", "Test $file file exists" );
    ok( -s $temp_dir."/edit_dir/$file", "Created new $file $file" );
    ok( File::Compare::compare( $data_dir."/edit_dir/$file", $temp_dir."/edit_dir/$file" ) == 0, "$file files match" );
}

done_testing();
