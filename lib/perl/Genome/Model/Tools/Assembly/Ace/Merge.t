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

my @test_aces = qw/ test_asm1.ace test_asm2.ace test_asm3.ace /;

foreach (@test_aces) {
    ok(-s $data_dir."/$_", "Test ace file $_ exists");
    symlink($data_dir."/$_", $temp_dir."/$_");
    ok(-s $temp_dir."/$_", "Linked $_");
}

my $temp_ace_list = $temp_dir.'/ace_list';
my $ace_list_fh = IO::File->new(">".$temp_ace_list) ||
    die "Can not create file handle to write ace list\n";
$ace_list_fh->print(map {$_."\n"} @test_aces);
$ace_list_fh->close;

ok(-s $temp_dir.'/ace_list', "Created temp ace list");
my $cmd = "gmt assembly ace merge --ace-list $temp_ace_list --directory $temp_dir";
ok(! system("$cmd"),"Ran successfully command:\n\ $cmd");

my $test_ace = $data_dir.'/edit_dir/merged.final.ace';
my $temp_ace = $temp_dir.'/merged.final.ace';

ok(-s $temp_ace, "Created merged ace file");
ok(-s $test_ace, "Test merged ace fiel exists");

my @diff = `sdiff -s $temp_ace $test_ace`;
is (scalar @diff, 0, "New merged ace file matches test merged ace file");

done_testing();

exit;
