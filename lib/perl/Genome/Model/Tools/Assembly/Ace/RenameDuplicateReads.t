#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::Assembly::Ace::RenameDuplicateReads' ) or die;

#check test suite dir
my $version = 'v1';
my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Assembly/RenameDuplicateReads-'.$version;
ok ( -d $test_dir, 'test suite dir' );

#make temp test dir
my $temp_dir = Genome::Sys->create_temp_directory();
ok( -d $temp_dir, 'made temp test dir' );
Genome::Sys->create_directory( $temp_dir.'/edit_dir' );
ok( -d $temp_dir.'/edit_dir', 'made temp edit_dir' );
Genome::Sys->create_directory( $temp_dir.'/phdball_dir' );
ok( -d $temp_dir.'/phdball_dir', 'made temp phdball dir' );

#copy test files
my $ace_file = 'test.ace.broke';
ok( -s $test_dir."/edit_dir/$ace_file", 'ace file exists in test dir' );
ok( File::Copy::copy( $test_dir."/edit_dir/$ace_file", $temp_dir."/edit_dir/$ace_file"), 'copied ace file' );
my $phdball = 'phdball.1';
ok( -s $test_dir."/phdball_dir/$phdball", 'phdball file exists in test dir' );
ok( File::Copy::copy( $test_dir."/phdball_dir/$phdball", $temp_dir."/phdball_dir/$phdball" ), 'copied phdball file' );

#create/execute tool
my $tool = Genome::Model::Tools::Assembly::Ace::RenameDuplicateReads->create(
    ace_in => $temp_dir."/edit_dir/$ace_file",
);
ok( $tool, 'created tool' );
ok( $tool->execute, 'executed tool' );

#compare output files
my $ace_out = 'test.ace.broke.duplicate_reads_renamed';
ok( -s $test_dir."/edit_dir/$ace_out", 'result ace file exists in test dir' );
ok( -s $temp_dir."/edit_dir/$ace_out", 'created new ace in temp dir' );
ok( File::Compare::compare( $test_dir."/edit_dir/$ace_out",$temp_dir."/edit_dir/$ace_out",) == 0, 'Result ace files matches' );

#<STDIN>;

done_testing();
