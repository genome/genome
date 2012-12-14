#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::Assembly::CreateOutputFiles::ContigsFromAce' ) or die;

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Assembly-CreateOutputFiles';
ok(-d $data_dir, "Found data directory: $data_dir") or die;

my $test_ace = $data_dir.'/edit_dir/velvet_asm.ace';
ok(-s $test_ace, "Test ace file exists") or die;

my $temp_dir = Genome::Sys->create_temp_directory();
Genome::Sys->create_directory( $temp_dir.'/edit_dir' );

ok(File::Copy::copy($test_ace, $temp_dir.'/edit_dir'),"Copied input ace file to temp dir") or die;

my $ec = system("chdir $temp_dir; gmt assembly create-output-files contigs-from-ace --acefile $temp_dir/edit_dir/velvet_asm.ace --directory $temp_dir");
ok($ec == 0, "Command ran successfully") or die;

foreach my $file ('contigs.bases', 'contigs.quals') {
    ok( -s $data_dir."/edit_dir/$file", "Test file $file exists");
    ok( -s $temp_dir."/edit_dir/$file", "New file $file created");
    ok( File::Compare::compare( $data_dir."/edit_dir/$file", $temp_dir."/edit_dir/$file" ) == 0, "$file files match");
}

done_testing();
