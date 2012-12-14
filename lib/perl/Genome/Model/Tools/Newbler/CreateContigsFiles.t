#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use above "Genome";
use File::Basename;
require File::Compare;

use_ok( 'Genome::Model::Tools::Newbler::CreateContigsFiles' ) or die;

#test suite dir
my $version = 'v1';
my $test_suite = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Newbler/CreateContigsFiles-'.$version;
ok( -d $test_suite, "Test suite dir exists" ) or die;

#temp test dir
my $temp_dir = Genome::Sys->create_temp_directory();
ok ( -d $temp_dir, "Temp test dir created" ) or die;

for my $test_type (qw/ scaffolded unscaffolded / ) {
    #check test suite for test type
    my $test_suite_dir = $test_suite."/$test_type";
    ok( -d $test_suite_dir, "Test suite dir exists for $test_type" ) or die;

    #temp test dir for test type
    my $temp_test_dir = $temp_dir."/$test_type";
    Genome::Sys->create_directory( $temp_test_dir );
    ok ( -d $temp_test_dir, "Sub dir for test type $test_type created" ) or die;
    
    #input files for each test type
    my @input_files;
    if ( $test_type eq 'scaffolded' ) {
        @input_files = &scaffolded_input_files;
    }
    elsif ( $test_type eq 'unscaffolded' ) {
        @input_files = &unscaffolded_input_files;
    }

    #link input files;
    for my $file ( @input_files ) {
        ok( -s $test_suite_dir."/$file", "Test suite $file file exists" ) or die;
        symlink( $test_suite_dir."/$file", $temp_test_dir."/$file" );
        ok( -l $temp_test_dir."/$file", "Linked $file file" ) or die;
    }
    
    #create/execute tool
    my $create = Genome::Model::Tools::Newbler::CreateContigsFiles->create(
        assembly_directory => $temp_test_dir,
        min_contig_length => 5000,
        );
    ok ($create, "Successfully created tool" ) or die;
    ok ($create->execute, "Successfully executed tool" ) or die;
    
    #compare output files
    for my $file ( 'contigs.bases', 'contigs.quals' ) {
        ok( -s $test_suite_dir."/consed/edit_dir/$file", "Test suite $file file exists" );
        ok( -s $temp_test_dir."/consed/edit_dir/$file", "Create $file file" );
        ok( File::Compare::compare( $test_suite_dir."/consed/edit_dir/$file", $temp_test_dir."/consed/edit_dir/$file" ) == 0, "$test_type $file files match" );
    }
    
}

#<STDIN>;

done_testing();



sub scaffolded_input_files {
    return (qw/ 454AllContigs.fna 454AllContigs.qual 454Scaffolds.txt / );
}

sub unscaffolded_input_files {
    return (qw/ 454AllContigs.fna 454AllContigs.qual / );
}
