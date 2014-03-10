#!/usr/bin/env perl

use above 'Genome';
use Data::Dumper;
use IO::File;
use Test::More;
use File::Basename qw/dirname/;
use File::Temp qw/tempdir/;

use strict;
use warnings;

my $pkg = 'Genome::Model::Tools::EdgeR::FilterCountsFile';

use_ok($pkg);

my $tmpdir = tempdir(CLEANUP => 1);
my $input_file = "$tmpdir/in.txt";
my $output_file = "$tmpdir/out.txt";
my $num_genes = 20;
my $num_normal = 2;
my $num_tumor = 2;

Genome::Model::Tools::EdgeR::Base::generate_counts_file($input_file, $num_genes, $num_normal, $num_tumor);

subtest "execute" => sub {
    my $cmd = $pkg->execute(
        counts_file             => $input_file,
        counts_per_million      => 60000,
        num_samples             => 4,
        filtered_counts_file    => $output_file,
    );

    my $in = Genome::Sys->open_file_for_reading($output_file);
    my $header = $in->getline;
    chomp $header;
    my @fields = split("\t", $header);
    is_deeply(["", "Normal1", "Normal2", "Tumor1", "Tumor2"], \@fields, "header is as expected");

    my $genome = $in->getline;
    ok($genome =~ m/^GENE0/, "GENE0 has passed filtering");

    my $next_line = $in->getline;
    ok(!defined($next_line), "GENE0 is the only gene that passed filtering");
};

done_testing();
