#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

use_ok ('Genome::Model::Tools::Soap::CreateSupercontigsAgpFile') or die;

#check test data files
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Soap/CreateSupercontigsAgpFile-v1';
ok(-d $data_dir, "Data dir exists") or die;

ok(-s $data_dir."/TEST.scafSeq", "Data dir test file exists") or die;

ok(-s $data_dir."/supercontigs.agp", "Data dir supercontigs.agp file exists") or die;

#create tmp test data dir
my $temp_dir = Genome::Sys->create_temp_directory();
ok(-d $temp_dir, "Temp test dir created") or die;

symlink( $data_dir.'/TEST.scafSeq', $temp_dir.'/TEST.scafSeq' );
ok( -s $temp_dir.'/TEST.scafSeq', "Linked TEST.scafSeq file" );

#create, execute command
my $create = Genome::Model::Tools::Soap::CreateSupercontigsAgpFile->create(
    assembly_directory => $temp_dir,
    min_contig_length => 50,
    );
ok($create, "Created gmt soap create-supercontigs-agp-file") or die;
ok(($create->execute) == 1, "Executed command") or die;
ok(-s $temp_dir.'/edit_dir/supercontigs.agp', "Created supercontigs.agp") or die;

#compare output files
ok(File::Compare::compare($data_dir.'/supercontigs.agp', $temp_dir.'/edit_dir/supercontigs.agp') ==0, "Output files match");

#print $temp_dir."\n"; <STDIN>;

done_testing();

exit;
