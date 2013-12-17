#!/usr/bin/env perl

use above 'Genome';
use Test::More;
use Data::Dumper;
use File::Temp qw/tempdir tempfile/;
use File::Slurp qw/read_file write_file/;

use strict;
use warnings;

my $pkg = 'Genome::Model::Tools::EdgeR::FilterStructures';
use_ok($pkg);

my $tempdir = tempdir(CLEANUP => 1);

my $counts = [
    [200, 200, 200, 200],
    [200, 200, 200, 100],
    [200, 200, 100, 100],
    [200, 100, 100, 100],
    [100, 100, 100, 100],
];
my $input_data = join("\n",
        join("\t", "GENE", map {"SAMPLE$_"} 0..3),
        make_matrix_text()
        );

my $counts_file = "$tempdir/counts.txt";

write_file($counts_file, $input_data);

subtest "Filter 25% >= 150" => make_test(25, 150, [0, 1]);
subtest "Filter 25% >= 200" => make_test(25, 200, [0]);
subtest "Filter 50% >= 150" => make_test(50, 150, [0, 1, 2]);

done_testing();

sub make_test {
    my ($percentile, $min_count, $expected_rows) = @_;
    return sub {
        my ($fh, $output_file) = tempfile(DIR => $tempdir);
        my $cmd = $pkg->create(
            counts_file => $counts_file,
            output_file => $output_file,
            percentile => $percentile,
            min_count => $min_count
            );
        ok($cmd, "created command");

        ok($cmd->execute, "executed command");
        my $expected_text = make_matrix_text(@$expected_rows);

        my @output = parse_output($output_file);
        my @expected = @{$counts}[@$expected_rows];
        is_deeply(\@output, \@expected, "Output is as expected")
            or diag("expected: " . Dumper(\@expected)
                . "\nactual: " . Dumper(\@output));
    }
}

sub make_matrix_text {
    my @indices = @_;
    if (!scalar @indices) {
        @indices = 0..$#$counts;
    }

    return join("\n",
        map {join("\t", "GENE$_", @{$counts->[$_]})} @indices);
}

sub parse_output {
    my $path = shift;
    my $output = read_file($path);

    my @lines = split("\n", $output);
    chomp @lines;
    shift @lines; # discard header

    return map {my @x = split("\t", $_); shift @x; \@x} @lines;
}
