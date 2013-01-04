#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::Velvet::ToConsed' ) or die;

my $test_suite_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Velvet-ToConsed';
ok ( -d $test_suite_dir, "Test suite exists" ) or die;

my $temp_test_dir = Genome::Sys->create_temp_directory();
ok ( -d $temp_test_dir, "Temp test dir created" ) or die;

#check necessary input files in test suite
for my $file (qw/ Sequences contigs.fa input.fastq velvet_asm.afg / ) {
    ok ( -s $test_suite_dir."/$file", "Test suite $file exists" );
    ok ( File::Copy::copy( $test_suite_dir."/$file", $temp_test_dir."/$file" ), "Copied $file to temp test dir" );
}

#create execute tool
my $cmd = Genome::Model::Tools::Velvet::ToConsed->create(
    assembly_directory => $temp_test_dir,
    no_phdball => 1, #otherwise this takes a looog time
    time_stamp => 'Mon Jun 20 15:11:40 2011',
);

ok ( $cmd, "Successfully created tool" ) or die;
ok ( $cmd->execute, "Successfully executed tool" ) or die;

#compare output file
my $test_suite_ace = $test_suite_dir.'/edit_dir/velvet_asm.ace';
my $temp_test_ace = $temp_test_dir.'/edit_dir/velvet_asm.ace';
ok ( -s $test_suite_ace, "Test suite ace file exists" );
ok ( -s $temp_test_ace, "Temp test ace file exists" );
ok ( File::Compare::compare( $test_suite_ace, $temp_test_ace ) == 1, "Ace files match except for one time stamp line" );

#check that phdball file did not get created since using --no-phdball option
ok ( ! -s $temp_test_dir.'/edit_dir/phd.ball.1', "Did not create phdball file" );

#<STDIN>;

done_testing();
