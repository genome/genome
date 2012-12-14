#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

my $module = 'Genome-Model-Tools-Assembly-Ace';
my $data_dir = $ENV{GENOME_TEST_INPUTS} . "/$module";

ok(-d $data_dir, "Found test data dir");
#create test temp dir
my $temp_dir = Genome::Sys->create_temp_directory();
ok(-d $temp_dir, "Made temp directory at $temp_dir");

#copy/link input files
 #ace files
my @test_aces = qw/ test_asm1.ace test_asm2.ace test_asm3.ace /;
foreach (@test_aces) {
    ok(-s $data_dir."/$_", "Test ace file $_ exists");
    symlink($data_dir."/$_", $temp_dir."/$_");
    ok(-s $temp_dir."/$_", "Linked $_");
}
 #ace list
my $data_ace_list = $data_dir.'/ace_list';
my $temp_ace_list = $temp_dir.'/ace_list';
ok(-s $data_ace_list, "Data dir ace_list exists");
symlink($data_ace_list, $temp_ace_list);
ok(-s $temp_ace_list, "Temp dir ace_list linked");
 #contigs list
my $data_contigs_list = $data_dir.'/export_contigs_list';
my $temp_contigs_list = $temp_dir.'/export_contigs_list';
ok(-s $data_contigs_list, "Data dir export contigs list exists");
symlink($data_contigs_list, $temp_contigs_list);
ok(-s $temp_contigs_list, "Temp dir export_contigs list linked");
#run test commands
my $cmd = "gmt assembly ace remove-contigs --ace-list $temp_ace_list --contigs-list $temp_contigs_list --directory $temp_dir";
ok(! system("$cmd"),"Ran successfully command:\n\ $cmd");

foreach (qw/ test_asm1.ace.contigs_removed.ace test_asm2.ace.contigs_removed.ace test_asm3.ace.contigs_removed.ace/) {
    ok(-s $data_dir."/$_", "Test $_ file exists");
    ok(-s $temp_dir."/$_", "$_ file get created");
    my @diff = `sdiff -s $data_dir/$_ $temp_dir/$_`;
    is(scalar @diff, 0, "$_ files match");
}

done_testing();
