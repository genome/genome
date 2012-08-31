#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use File::Temp;

my $abyss_version = '1.2.7';
my $class = "Genome::Model::Tools::Abyss::Parallel";
use_ok($class);

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my %create_params = (
    version => $abyss_version,
    kmer_size => 50,
    min_pairs => 11,
    num_jobs => 16,
    name => 'test',
    fastq_a => "$tmpdir/a.fq",
    fastq_b => "$tmpdir/b.fq",
    output_directory => $tmpdir,
);

# valid kmer_sizes
my @rv = $class->get_kmer_sizes(7);
is_deeply(\@rv, [7], 'simple integer works');

@rv = $class->get_kmer_sizes("1..3"); 
is_deeply(\@rv, [1,2,3], 'range works');

@rv = $class->get_kmer_sizes("10..20:2");
is_deeply(\@rv, [10,12,14,16,18,20], 'range with step works');

@rv = $class->get_kmer_sizes("4,5,6");
is_deeply(\@rv, [4,5,6], 'list works');

@rv = $class->get_kmer_sizes("1,2,3,10..13,20");
is_deeply(\@rv, [1,2,3,10,11,12,13,20], 'list with range works');

@rv = $class->get_kmer_sizes("1,2,3,10..20:3,99");
is_deeply(\@rv, [1,2,3,10,13,16,19,99], 'list with range+step works');

@rv = $class->get_kmer_sizes("25,40..50:2");
is_deeply(\@rv, [25,40,42,44,46,48,50], 'list with range+step works');


# invalid kmer_sizes
eval { $class->get_kmer_sizes("cat"); };
ok($@, "cat is not a valid number");

eval { $class->get_kmer_sizes("9.3"); };
ok($@, "floating point numbers not accepted");

eval { $class->get_kmer_sizes("9..3"); };
ok($@, "range where end < start is error");

eval { $class->get_kmer_sizes("3..9:0"); };
ok($@, "range where step = 0 is error");

eval { $class->get_kmer_sizes("3..9:-1"); };
ok($@, "range where step < 0 is error");

# object creation
eval { $class->create(); };
ok($@, 'create with no params fails');

my $obj = $class->create(%create_params);
ok($obj, 'created object');

done_testing();
