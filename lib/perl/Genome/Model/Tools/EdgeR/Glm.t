#!/usr/bin/env perl

use above 'Genome';
use Data::Dumper;
use IO::File;
use Test::More;
use File::Basename qw/dirname/;
use File::Temp qw/tempdir/;

use strict;
use warnings;

my $pkg = 'Genome::Model::Tools::EdgeR::Glm';

use_ok($pkg);

my $tmpdir = tempdir(CLEANUP => 1);
my $input_file = "$tmpdir/in.txt";
my $output_file = "$tmpdir/out.txt";
my $num_genes = 20;
my $num_normal = 2;
my $num_tumor = 2;

Genome::Model::Tools::EdgeR::Base::generate_counts_file($input_file, $num_genes, $num_normal, $num_tumor);

subtest "subject only used once" => sub {
        my $cmd = $pkg->create(
        counts_file => $input_file,
        groups => "A,A,B,B",
        subjects => "s1,s2,s3,s4",
        output_file => "/dev/null"
    );

    my $rv = 0;
    eval {
        $rv = $cmd->_validate_params();
    };
    ok($@, "Subject only used once is an error");
    ok(!$rv);
};

subtest "only one subject" => sub {
        my $cmd = $pkg->create(
        counts_file => $input_file,
        groups => "A,A,B,B",
        subjects => "s1,s1,s1,s1",
        output_file => "/dev/null"
    );

    my $rv = 0;
    eval {
        $rv = $cmd->_validate_params();
    };
    ok($@, "Only one subject is an error");
    ok(!$rv);
};

subtest "different number of groups and subjects" => sub {
        my $cmd = $pkg->create(
        counts_file => $input_file,
        groups => "A,A,B,B",
        subjects => "s1,s1,s1",
        output_file => "/dev/null"
    );

    my $rv = 0;
    eval {
        $rv = $cmd->_validate_params();
    };
    ok($@, "Different nuber of groups and subjects is an error");
    ok(!$rv);
};

subtest "each subject is not in each group" => sub {
        my $cmd = $pkg->create(
        counts_file => $input_file,
        groups => "A,A,B,B",
        subjects => "s1,s1,s2,s2",
        output_file => "/dev/null"
    );

    my $rv = 0;
    eval {
        $rv = $cmd->_validate_params();
    };
    ok($@, "Each subject is not in each group is an error");
    ok(!$rv);
};

subtest "execute" => sub {
    my $cmd = $pkg->create(
        counts_file => $input_file,
        groups => 'normal,normal,tumor,tumor',
        subjects => 's1,s2,s1,s2',
        output_file => $output_file,
    );
    ok($cmd, "Created command");
    ok($cmd->execute, "Executed command");

    my $in = Genome::Sys->open_file_for_reading($output_file);
    my $header = $in->getline;
    chomp $header;
    my @fields = split("\t", $header);
    is_deeply(["", "genes", "logFC", "logCPM", "LR", "PValue", "FDR", "test.result"], \@fields, "header is as expected");
    my %gene_results;
    my %gene_tests;
    while (my $line = $in->getline) {
        chomp $line;
        my ($coeff, $name, $x, $y, $z, $pvalue, $adjusted_pvalue, $test) = split("\t", $line);
        $gene_results{$name}{pvalue} = $adjusted_pvalue;
        $gene_results{$name}{test} = $test;
    }

    my @missing = grep {
            my $n = "GENE$_";
            !exists $gene_results{$n}{pvalue} or $gene_results{$n}{pvalue} eq 'NA'
        } 0..$num_genes - 1;

    is(scalar @missing, 0, "Got data for all genes") or diag("Missing: " . Dumper(\@missing));

    ok($gene_results{GENE0}{pvalue} < 0.001,
        "GENE0 is differentially expressed with high confidence");

    is($gene_results{GENE0}{test}, 1,
        "GENE0 passed significance test for DE");

    my @fail = grep { $gene_results{"GENE$_"}{pvalue} < 0.5 } 1..$num_genes - 1;
    is(scalar @fail, 0, "All other genes are p >= 0.5") or diag("Failed: " . Dumper(\@fail));

    @fail = grep { $gene_results{"GENE$_"}{test} != 0 } 1..$num_genes - 1;
    is(scalar @fail, 0, "All other genes fail significance test") or diag("Failed: " . Dumper(\@fail));

};

done_testing();



