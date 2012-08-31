#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use above "Genome";
require File::Compare;

my $archos = `uname -a`;
unless ($archos =~ /64/) {
    plan skip_all => "Must run from 64-bit machine";
}

use_ok( 'Genome::Model::Tools::Newbler::ToPcapAce' ) or die;

#test suite dir
my $version = 'v2';
my $test_suite = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Newbler/ToPcapAce-'.$version;
ok( -d $test_suite, "Test suite dir exists" ) or die;

my $temp_dir = Genome::Sys->create_temp_directory();
ok( -d $temp_dir, "Temp test dir creates" );

for my $test_type ( qw/ scaffolded unscaffolded / ) {
    #test suite type dir
    my $suite_type_dir = $test_suite."/$test_type";
    ok( -d $suite_type_dir, "Test suite dir exists for $test_type test" );

    #temp test dir for type
    my $temp_type_dir = $temp_dir."/$test_type";
    Genome::Sys->create_directory( $temp_type_dir );
    ok( -d $temp_type_dir, "Created temp dir for $test_type test" );
    Genome::Sys->create_directory( $temp_type_dir.'/consed' );
    ok( -d $temp_type_dir.'/consed', "Made consed dir" );
    Genome::Sys->create_directory( $temp_type_dir.'/consed/edit_dir' );
    ok( -d $temp_type_dir.'/consed/edit_dir', "Made edit_dir" );

    #input files for each test
    my @input_files;
    if ( $test_type eq 'scaffolded' ) {
        @input_files = &scaffolded_input_files;
    } else {
        @input_files = &unscaffolded_input_files;
    }
    ok( @input_files == 2, "Got two input files" );

    #link input files
    for my $file ( @input_files ) {
        ok( -s $suite_type_dir."/$file", "Test suite $file file exists" );
        symlink( $suite_type_dir."/$file", $temp_type_dir."/$file" );
        ok( -l $temp_type_dir."/$file", "Linked $file file" );
    }

    #create/execute tool
    my $create = Genome::Model::Tools::Newbler::ToPcapAce->create(
        assembly_directory => $temp_type_dir,
        min_contig_length => 200,
    );
    ok( $create, "Created tool" );
    ok( $create->execute, "Executed tool" );

    #compare output files
    for my $file ('/consed/edit_dir/Pcap.454Contigs.ace', '/consed/edit_dir/gap.txt' ) {
        ok( -s $suite_type_dir."/$file", "Test suite $file file exists" );
        ok( -s $temp_type_dir."/$file", "Test created $file file" );
        ok( File::Compare::compare( $suite_type_dir."/$file",$temp_type_dir."/$file") == 0, "$file files match" );
    }


}

#<STDIN>;

done_testing();

exit;

sub scaffolded_input_files {
    return ( '454Scaffolds.txt', '/consed/edit_dir/454Contigs.ace.1' );
}

sub unscaffolded_input_files {
    return ( '454AllContigs.fna', '/consed/edit_dir/454Contigs.ace.1' );
}
