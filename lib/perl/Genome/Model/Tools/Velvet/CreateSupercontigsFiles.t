#!/usr/bin/env genome-perl

use strict;
use warnings;

require File::Compare;
use above "Genome";
use Test::More;

use_ok ( 'Genome::Model::Tools::Velvet::CreateSupercontigsFiles' ) or die;

my $version = 'v3';
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Velvet/CreateSupercontigsFiles-'.$version;
ok(-d $data_dir, "Found data directory: $data_dir") or die;

my $temp_dir = Genome::Sys->create_temp_directory();
ok( -d $temp_dir, "Temp test dir created" ) or die;

#link project dir files
my $file = 'velvet_asm.afg';
ok(-s $data_dir."/$file", "Data dir $file file exists"); 
symlink ($data_dir."/$file", $temp_dir."/$file");
ok(-s $temp_dir."/$file", "Tmp dir file file exists");

#create / execute tool
my $create = Genome::Model::Tools::Velvet::CreateSupercontigsFiles->create (
    assembly_directory => $temp_dir,
    min_contig_length => 50,
    );
ok( $create, "Created tool");
ok( $create->execute, "Successfully executed tool");

foreach ('supercontigs.fasta', 'supercontigs.agp') {
    ok(-s $data_dir."/edit_dir/$_", "Data dir $_ file exists");
    ok(-s $temp_dir."/edit_dir/$_", "Tmp dir $_ file got created");
    ok(File::Compare::compare($data_dir."/edit_dir/$_", $temp_dir."/edit_dir/$_") == 0, "$_ files match");
}

#<STDIN>;

done_testing();
