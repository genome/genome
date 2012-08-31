#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::Assembly::CreateOutputFiles::ReadInfo' ) or die;

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Assembly-CreateOutputFiles';
ok(-d $data_dir, "Found data directory: $data_dir") or die;

my $test_ace = $data_dir.'/edit_dir/velvet_asm.ace';
ok(-s $test_ace, "Found test ace file") or die;

#make test dir
my $temp_dir = Genome::Sys->create_temp_directory();
Genome::Sys->create_directory( $temp_dir.'/edit_dir' );

ok(File::Copy::copy($test_ace, $temp_dir.'/edit_dir'),"Copied input ace file to temp dir");

my $ec = system("chdir $temp_dir; gmt assembly create-output-files read-info --acefile $temp_dir/edit_dir/velvet_asm.ace --directory $temp_dir");
ok($ec == 0, "Command ran successfully");

ok( -s $data_dir.'/edit_dir/readinfo.txt', "Test readinfo file exists" );
ok( -s $temp_dir.'/edit_dir/readinfo.txt', "New readinfo file created" );
ok( File::Compare::compare( $data_dir.'/edit_dir/readinfo.txt', $temp_dir.'/edit_dir/readinfo.txt' ) == 0, "Files match" );

done_testing();

exit;
