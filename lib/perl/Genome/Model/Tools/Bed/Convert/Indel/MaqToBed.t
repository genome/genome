#!/usr/bin/env genome-perl

use strict;
use warnings;


use Test::More tests => 5;

use above 'Genome';

use_ok('Genome::Model::Tools::Bed::Convert::Indel::MaqToBed');

my $tmpdir = File::Temp::tempdir('Bed-Convert-Indel-MaqToBedXXXXX', CLEANUP => 1, TMPDIR => 1);
my $output_file = join('/', $tmpdir, 'output');

my $input_file = __FILE__ . '.input';
my $expected_file = __FILE__ . '.expected';

my $command = Genome::Model::Tools::Bed::Convert::Indel::MaqToBed->create( source => $input_file, output => $output_file );
ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(!$command->_need_chrom_sort, 'ordered chromosomes in input, no sort required');

my $diff = Genome::Sys->diff_file_vs_file($output_file, $expected_file);
ok(!$diff, 'output matched expected result')
    or diag("diff results:\n" . $diff) && diag('Note that if you need to update the expected file you may also need to replace the expected directory for the long test of `gmt detect-variants maq`');
