#!/usr/bin/env genome-perl

use strict;
use warnings;

require File::Compare;
use above "Genome";
use Test::More;

use_ok( 'Genome::Model::Tools::Velvet::CreateContigsFiles' );

#TODO - move to correct test suite module dir when all tests are configured
my $version = 'v3';
my $data_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Velvet/CreateReadsFiles-".$version;

ok(-d $data_dir, "Found data directory: $data_dir") or die;

my $temp_dir = Genome::Sys->create_temp_directory();

#make edit_dir
mkdir $temp_dir.'/edit_dir';
ok(-d $temp_dir.'/edit_dir', "Made edit_dir in temp dir");

#link project dir files
foreach ('velvet_asm.afg', 'Sequences', 'contigs.fa') {
    ok(-s $data_dir.'/'.$_, "Data dir $_ file exists"); 
    symlink ($data_dir.'/'.$_, $temp_dir.'/'.$_);
    ok(-s $temp_dir.'/'.$_, "Tmp dir $_ file exists");
}

#link edit_dir files
ok(-s $data_dir.'/edit_dir/gap.txt', "Data dir gap.txt file exists");
symlink($data_dir.'/edit_dir/gap.txt', $temp_dir.'/edit_dir/gap.txt');
ok(-s $temp_dir.'/edit_dir/gap.txt', "Tmp edit_dir gap.txt file exists");

my $create = Genome::Model::Tools::Velvet::CreateReadsFiles->create(
    assembly_directory => $temp_dir,
    min_contig_length => 50,
    );
ok( $create, "Created tool");
ok( $create->execute, "Successfully executed tool" );

foreach ('readinfo.txt', 'reads.placed') {
    ok(-s $data_dir."/edit_dir/$_", "Data dir $_ file exists");
    ok(-s $temp_dir."/edit_dir/$_", "Tmp dir $_ file got created");
    ok(File::Compare::compare($data_dir."/edit_dir/$_", $temp_dir."/edit_dir/$_") == 0, "$_ files match");
}

#<STDIN>;

done_testing();

exit;
