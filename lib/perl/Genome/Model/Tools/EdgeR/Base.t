#!/usr/bin/env perl

use above 'Genome';
use Data::Dumper;
use IO::File;
use Test::More;
use File::Basename qw/dirname/;
use File::Temp qw/tempdir/;

use strict;
use warnings;

my $pkg = 'Genome::Model::Tools::EdgeR::Base';

use_ok($pkg);

my $tmpdir = tempdir(CLEANUP => 1);
my $input_file = "$tmpdir/in.txt";
my $output_file = "$tmpdir/out.txt";
my $num_genes = 20;
my $num_normal = 1;
my $num_tumor = 2;

Genome::Model::Tools::EdgeR::Base::generate_counts_file($input_file, $num_genes, $num_normal, $num_tumor);

subtest "bad p-values" => sub {
    my @bad_pvalues = (-1, 0, 1, 2);

    for my $p (@bad_pvalues) {
        my $cmd = $pkg->create(
            counts_file => $input_file,
            groups => "1,2,2",
            output_file => "/dev/null",
            p_value => $p,
        );

        my $rv = 0;
        eval {
            $rv = $cmd->_validate_params();
        };
        ok($@, "p_value $p is invalid");
        ok(!$rv);
    }
};

subtest "no replication" => sub {
    my $cmd = $pkg->create(
        counts_file => $input_file,
        groups => "1,2,3",
        output_file => "/dev/null"
    );

    my $rv = 0;
    eval {
        $rv = $cmd->_validate_params();
    };
    ok($@, "No replication is an error");
    ok(!$rv);
};

done_testing();
