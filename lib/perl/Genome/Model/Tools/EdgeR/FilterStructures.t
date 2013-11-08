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

my $counts = [
    [200, 200, 200, 200],
    [200, 200, 200, 100],
    [200, 200, 100, 100],
    [200, 100, 100, 100],
    [100, 100, 100, 100],
];

my $data = join("\n",
    join("\t", "GENE", map {"SAMPLE$_"} 0..3),
    map {join("\t", "GENE$_", @{$counts->[$_]})} 0..$#$counts);

my $tempdir = tempdir(CLEANUP => 1);
my $counts_file = "$tempdir/counts.txt";

write_file($counts_file, $data);

subtest "Filter 75% >= 150" => make_test(75, 150, ["GENE0", "GENE1"]);
subtest "Filter 75% >= 200" => make_test(75, 200, ["GENE0"]);
subtest "Filter 50% >= 150" => make_test(50, 150, ["GENE0", "GENE1", "GENE2"]);

done_testing();

sub make_test {
    my ($percentile, $min_count, $expected) = @_;
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

        my @genes = get_genes($output_file);
        is_deeply($expected, \@genes, "Expected genes passed filter");
    }
}

sub get_genes {
    my $path = shift;
    my $output = read_file($path);

    my @lines = split("\n", $output);
    chomp @lines;
    shift @lines; # discard header

    return map {my @f = split("\t", $_); $f[0]} @lines;
}
