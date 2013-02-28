#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 9;

use above 'Genome';

use_ok('Genome::Model::Tools::Bed::Convert::Snv::SniperToBed');

my $tmpdir = File::Temp::tempdir('Bed-Convert-Snv-SniperToBedXXXXX', CLEANUP => 1, TMPDIR => 1);

my $input_file = __FILE__ . '.input';
my $expected_file = __FILE__ . '.expected';

verify_output($input_file, "$tmpdir/10col.out", $expected_file);
verify_output("$input_file-26", "$tmpdir/26col.out", "$expected_file-26");

sub verify_output {
    my ($input_file, $output_file, $expected_file) = @_;
    my $command = Genome::Model::Tools::Bed::Convert::Snv::SniperToBed->create( source => $input_file, output => $output_file );
    ok($command, 'Command created');
    my $rv = $command->execute;
    ok($rv, 'Command completed successfully');
    ok(-s $output_file, "output file created");

    my $diff = Genome::Sys->diff_file_vs_file($output_file, $expected_file);
    ok(!$diff, 'output matched expected result')
        or diag("diff results:\n" . $diff);
}
