#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Genome::Utility::Test;
use Test::More tests => 5;

my $input_path = File::Basename::dirname(__FILE__) . '/../Export/Metadata.t.expected-output';
ok(-e $input_path, "input path exists: $input_path");

my $expected_log_path = __FILE__ . '.expected-output';
ok(-e $expected_log_path, "expected log output file $expected_log_path exists");

my $actual_log_path = Genome::Sys->create_temp_file_path(); 
$ENV{UR_DBI_NO_COMMIT} = 1;

if ($ARGV[0] && $ARGV[0] eq 'REBUILD') {
    unlink $expected_log_path;
    $actual_log_path = $expected_log_path;
    warn "regenerating $expected_log_path...";
}

my $result = Genome::Model::Command::Import::Metadata->execute(input_path => $input_path, log_path => $actual_log_path, verbose => 1, skip_file_db_install => 1);
ok($result, "ran");
ok(-e $actual_log_path, "actual_log_path $actual_log_path exists");

Genome::Utility::Test::compare_ok(
    $actual_log_path, 
    $expected_log_path,
    'log matches',
    filters => [
        sub{ 
            my $line = shift;
            return if $line =~ /total_kb/i; # ignore disk volumes
            return $line;
        },
    ],
);

done_testing(5);
