#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Utility::Test;
use File::Spec;
use File::Temp;
use Test::More tests => 5;

my $pkg = 'Genome::Model::Tools::Sx::Split::ByNs';
use_ok($pkg);

my $dir = File::Spec->join( Genome::Config::get('test_inputs'), 'Genome-Model-Tools-Sx', 'SplitByNs');
my $in_fasta = File::Spec->join($dir, 'in.fasta');
ok(-s $in_fasta, 'in fasta');
my $expected_fasta = File::Spec->join($dir, 'expected.fasta');
ok(-s $expected_fasta, 'expected fasta');

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $out_fasta = File::Spec->join($tmp_dir, 'out.fasta');

my $splitter = $pkg->execute(
    input  => [ $in_fasta ],
    output => [ $out_fasta ],
    number_of_ns => 10,
);
ok($splitter->result, 'execute');
Genome::Utility::Test::compare_ok($out_fasta, $expected_fasta, 'output fasta matches expected');

done_testing();
