#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::Assembly::CreateOutputFiles::SupercontigsAgp' ) or die;

my $data_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Assembly/CreateSupercontigsAgp-v1";
ok(-d $data_dir, "Found data directory: $data_dir") or die;

#make test dir
my $temp_dir = Genome::Sys->create_temp_directory();
Genome::Sys->create_directory( $temp_dir.'/edit_dir' );

#copy files
for my $file ('gap.txt', 'contigs.bases') {
    ok( -e $data_dir.'/edit_dir/'.$file, "Test $file exists" );
    ok( File::Copy::copy($data_dir."/edit_dir/$file", $temp_dir."/edit_dir/$file"), "Copied $file to temp dir" );
}

#run
my $ec = system("chdir $temp_dir; gmt assembly create-output-files supercontigs-agp --directory $temp_dir");
ok($ec == 0, "Command ran successfully") or die;

#compare files
ok( -s $data_dir.'/edit_dir/supercontigs.agp', "Data dir supercontigs.agp file exists" );
ok( -s $temp_dir.'/edit_dir/supercontigs.agp', "Supercontigs.agp file created" );
ok( File::Compare::compare( $temp_dir.'/edit_dir/supercontigs.agp', $data_dir.'/edit_dir/supercontigs.agp' ) == 0, "Files match" );

#print $temp_dir."\n"; <STDIN>;
done_testing();

exit;

