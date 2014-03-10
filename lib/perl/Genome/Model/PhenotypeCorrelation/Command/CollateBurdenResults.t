#!/usr/bin/env perl

use above 'Genome';
use JSON;
use Test::More;
use File::Temp qw(tempdir);
use File::Slurp qw(write_file read_file);

use strict;
use warnings;

my $pkg = 'Genome::Model::PhenotypeCorrelation::Command::CollateBurdenResults';
use_ok($pkg);

my @test_names = qw(A B C);
my $expected_header = join(",", qw(Trait Gene), @test_names);

sub generate_data {
    my @genes = qw(G1 G2);
    my @phenos = qw(P1 P2);
    my @lines;
    my $counter = 0;
    my $header = "Trait,Gene,A,B,C";
    for my $g (@genes) {
        for my $p (@phenos) {
            my @tests = map {($counter + $_) / 10.0} 1..3;
            $counter += 3;
            push @lines, [$p, $g, @tests];
        }
    }
    return $header, @lines;
}

my ($header, @lines) = generate_data();

my $expected = join("\n", $header, sort map {join(",", @$_)} @lines) . "\n";

my $tmpdir = tempdir(CLEANUP => 1);

my $output_file = join("/", $tmpdir, "output.csv");
my $null_file = join("/", $tmpdir, "T3_Gene3.null");
write_file($null_file, "");
my %genes;
my %phenos;
for my $l (@lines) {
    $phenos{$l->[0]} = 1;
    $genes{$l->[1]} = 1;
    my $filename = sprintf("%s/%s_%s.burden.csv", $tmpdir, $l->[0], $l->[1]);
    write_file($filename, join("\n", $header, join(",", @$l)));
}

my @job_params;
for my $g (keys %genes) {
    for my $p (keys %phenos) {
        push(@job_params, {gene => $g, phenotype => $p});
    }
}

my $json = new JSON;
my @encoded_params = map {$json->encode($_)} @job_params;
my $cmd = $pkg->create(
    results_directory => $tmpdir,
    output_file => $output_file,
    job_params => \@encoded_params
    );

ok($cmd, "created command");
ok($cmd->execute, "executed command");
ok(-s $output_file, "output file exists and is nonempty");

my $data = read_file($output_file);
is($data, $expected, "got expected result");

done_testing();
