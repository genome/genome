#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 5;

my $input_path = File::Basename::dirname(__FILE__) . '/../Export/Metadata.t.expected-output';
ok(-e $input_path, "input path exists: $input_path");

my $expected_log_path = __FILE__ . '.expected-output';
ok(-e $expected_log_path, "expected log output file $expected_log_path exists");

my $actual_log_path = Genome::Sys->create_temp_file_path(); 
$ENV{UR_DBI_NO_COMMIT} = 1;

my $result = Genome::Model::Command::Import::Metadata->execute(input_path => $input_path, log_path => $actual_log_path);
ok($result, "ran");
ok(-e $actual_log_path, "actual_log_path $actual_log_path exists");

my @diff = grep { $_ !~ /total_kb/ } `sdiff -s -w 500 $expected_log_path $actual_log_path`;
is(scalar(@diff), 0, "no differences")
    or diag(@diff);



