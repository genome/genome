#!/usr/bin/env perl

use strict;
use warnings;

use Test::More tests => 13;

use above 'Genome';

use_ok('Genome::Model::Tools::Bed::Convert::Snv::StrelkaToBed');

my $tmpdir = File::Temp::tempdir('Bed-Convert-Snv-StrelkaToBedXXXXX', CLEANUP => 1, TMPDIR => 1);

my $input_file = __FILE__ . '.input';
my $expected_file = __FILE__ . '.expected';

my $output_file = "$tmpdir/actual.out";

# all

my $command = Genome::Model::Tools::Bed::Convert::Snv::StrelkaToBed->create( source => $input_file, output => $output_file );
ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file created");

my $diff = Genome::Sys->diff_file_vs_file($output_file, $expected_file);
ok(!$diff, 'output matched expected result')
    or diag("diff results:\n" . $diff);

# hq

$command = Genome::Model::Tools::Bed::Convert::Snv::StrelkaToBed->create( source => $input_file, output => $output_file . '.hqonly', limit_variants_to => 'hq' );
ok($command, 'Command created');
$rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file created");

$diff = Genome::Sys->diff_file_vs_file($output_file . '.hqonly', $expected_file . '.hqonly');
ok(!$diff, 'output matched expected result')
    or diag("diff results:\n" . $diff);

# lq

$command = Genome::Model::Tools::Bed::Convert::Snv::StrelkaToBed->create( source => $input_file, output => $output_file . '.lqonly', limit_variants_to => 'lq' );
ok($command, 'Command created');
$rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file created");

$diff = Genome::Sys->diff_file_vs_file($output_file . '.lqonly', $expected_file . '.lqonly');
ok(!$diff, 'output matched expected result')
    or diag("diff results:\n" . $diff);

#print $tmpdir, "\n";
