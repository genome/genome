#!/usr/bin/env genome-perl

use strict;
use warnings;


use Test::More tests => 4;

use above 'Genome';

use_ok('Genome::Model::Tools::Bed::Convert::Indel::VarscanSomaticToBed');

my $tmpdir = File::Temp::tempdir('Bed-Convert-Indel-VarscanSomaticToBedXXXXX', DIR => "$ENV{GENOME_TEST_TEMP}/", CLEANUP => 1);
my $output_file = join('/', $tmpdir, 'output');

my $input_file = __FILE__ . '.input';
my $expected_file = __FILE__ . '.expected';

# TODO FIX THIS to actually test VarscanSomatic instead of Varscan. Requires fixing up test data as well
my $command = Genome::Model::Tools::Bed::Convert::Indel::VarscanToBed->create( source => $input_file, output => $output_file );
ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');

my $diff = Genome::Sys->diff_file_vs_file($output_file, $expected_file);
ok(!$diff, 'output matched expected result')
    or diag("diff results:\n" . $diff);
