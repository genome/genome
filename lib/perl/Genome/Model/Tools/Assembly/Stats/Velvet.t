#!/usr/bin/env genome-perl

use strict;
use warnings;

use Cwd;
use above "Genome";
use Test::More;


my $data_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Assembly-Stats/Velvet_v3";
ok(-d $data_dir, "Found data directory: $data_dir");

#create temp test dir
my $temp_dir = Genome::Sys->create_temp_directory();

#make edit_dir in temp_dir
mkdir $temp_dir.'/edit_dir';
ok(-d $temp_dir.'/edit_dir', "Made edit_dir in temp test_dir");

#link assembly dir files needed to run stats
my @files_to_link = qw/ velvet_asm.afg Sequences /;
foreach (@files_to_link) {
    ok(-s $data_dir."/$_", "Test data file $_ exists");
    symlink($data_dir."/$_", $temp_dir."/$_");
    ok(-s $temp_dir."/$_", "Linked file $_ in tmp test dir");
}

#link assembly/edit_dir files needed to run stats
@files_to_link = qw/ velvet_asm.ace test.fasta.gz test.fasta.qual.gz
                        contigs.bases contigs.quals reads.placed readinfo.txt /;
foreach my $file (@files_to_link) {
    ok(-s $data_dir."/edit_dir/$file", "Test data file $file file exists");
    symlink($data_dir."/edit_dir/$file", $temp_dir."/edit_dir/$file");
    ok(-s $temp_dir."/edit_dir/$file", "Linked file $file in tmp test dir"); 
}

#create stats
ok(system ("gmt assembly stats velvet --assembly-directory $temp_dir --out-file $temp_dir/edit_dir/stats.txt --no-print-to-screen") == 0, "Command ran successfully");

#check for stats files
my $temp_stats = $temp_dir.'/edit_dir/stats.txt';
my $data_stats = $data_dir.'/edit_dir/stats.txt';

ok(-s $temp_stats, "Tmp test dir stats.txt file exists");
ok(-s $data_stats,, "Test data dir stats.txt file exists");

#compare files
my @diff = `sdiff -s $data_stats $temp_stats`;
is(scalar @diff, 0, "Stats files match") or diag(@diff);

done_testing();


 
