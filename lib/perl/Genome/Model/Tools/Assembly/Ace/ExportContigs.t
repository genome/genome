#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
require File::Compare;

my $module = 'Genome-Model-Tools-Assembly-Ace';
my $data_dir = $ENV{GENOME_TEST_INPUTS} . "/$module";

ok(-d $data_dir, "Found test data dir");

my $temp_dir = Genome::Sys->create_temp_directory();
ok(-d $temp_dir, "Made temp directory at $temp_dir");

my @test_aces = qw/ test_asm1.ace test_asm2.ace test_asm3.ace/;

foreach (@test_aces) {
    ok(-s $data_dir."/$_", "Test dir file $_ exists");
    symlink($data_dir."/$_", $temp_dir."/$_");
    ok(-s $temp_dir."/$_", "Linked $_");
}

my $ace_list_fh = IO::File->new(">".$temp_dir.'/ace_list') ||
    die "Can not create file handle to write ace list\n";
$ace_list_fh->print(map {$_."\n"} @test_aces);
$ace_list_fh->close;
ok(-s $temp_dir.'/ace_list', "Created temp ace list");

ok(-s $data_dir.'/export_contigs_list', "Test export contigs list exists");
symlink($data_dir.'/export_contigs_list', $temp_dir.'/export_contigs_list');
ok(-s $temp_dir.'/export_contigs_list', "Linked export contigs list to temp dir");

my $cmd = "gmt assembly ace export-contigs --ace-list $temp_dir/ace_list --contigs-list $temp_dir/export_contigs_list --directory $temp_dir";
ok (! system("$cmd"), "Ran command successfully");

my @result_aces = qw/ test_asm1.ace.exported_contigs.ace test_asm3.ace.exported_contigs.ace
                   test_asm2.ace.exported_contigs.ace /;
foreach (@result_aces) {
    ok(-s $data_dir."/$_", "Data $_ file exists");
    ok(-s $temp_dir."/$_", "Test created $_ file");
    my @diffs = `sdiff -s $data_dir/$_ $temp_dir/$_`;
    is (scalar @diffs, 0, "$_ files match");
}

done_testing();

exit;
