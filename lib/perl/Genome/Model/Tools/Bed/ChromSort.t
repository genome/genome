#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 8;
use above 'Genome';
use_ok('Genome::Model::Tools::Bed::ChromSort');

my $tmpdir = File::Temp::tempdir('Bed-ChromSort-XXXXX', CLEANUP => 1, TMPDIR => 1);
my $output_file = join('/', $tmpdir, 'output');
my $input_file = __FILE__ . '.input';
my $expected_file = __FILE__ . '.expected';

my $cmd = Genome::Model::Tools::Bed::ChromSort->create(
    input => $input_file,
    output => $output_file
    );

ok($cmd, 'Created command');
ok($cmd->execute(), 'Executed command');

my $diff = Genome::Sys->diff_file_vs_file($output_file, $expected_file);
ok(!$diff, 'output matched expected result') or diag("diff results:\n" . $diff);

#test the new functionality with skipping lines
my $input_with_header = $output_file . ".with_cat";
my $output_with_header = $input_with_header . ".after_sort";
ok(Genome::Sys->shellcmd(cmd => "echo 'header1\nheader2' | cat - $output_file >$input_with_header"));

my $cmd2 = Genome::Model::Tools::Bed::ChromSort->create(
    input => $input_with_header,
    output => $output_with_header,
    skip_lines => 2,
);

ok($cmd2, 'Created command for header test');
ok($cmd2->execute(), 'Executed command');

my $diff2 = Genome::Sys->diff_file_vs_file($input_with_header, $output_with_header);
ok(!$diff2, 'output with header matched expected result') or diag("diff results:\n" . $diff2);
   $DB::single = 1; 

