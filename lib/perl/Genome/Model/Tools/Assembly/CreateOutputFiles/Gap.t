#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;

use_ok( 'Genome::Model::Tools::Assembly::CreateOutputFiles::Gap' ) or die;
              
my $data_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Assembly-CreateOutputFiles";
ok(-d $data_dir, "Found data directory: $data_dir") or die;

#test gap.txt file
my $test_contigs_file = $data_dir.'/contigs.fa';
ok(-s $test_contigs_file, "Found test contigs.fa file");

#create temp test dir
my $temp_dir = Genome::Sys->create_temp_directory();
Genome::Sys->create_directory( $temp_dir.'/edit_dir' );

#copy input file
ok(File::Copy::copy($test_contigs_file, $temp_dir),"Copied input contigs file");

#run
my $ec = system("chdir $temp_dir; gmt assembly create-output-files gap --directory $temp_dir");
ok($ec == 0, "Command ran successfully");

#check, compare output files
ok( -e $data_dir.'/edit_dir/gap.txt', "Test gap.txt file exists" );
ok( -e $temp_dir.'/edit_dir/gap.txt', "New gap.txt file created" ); #can be zero size
my @diffs = `sdiff -s $data_dir/edit_dir/gap.txt $temp_dir/edit_dir/gap.txt`;
is(scalar (@diffs), 0, "New gap file matches test gap file");

done_testing();

exit;
