#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Sys;
use Genome::Utility::Test;
use File::Spec;
use File::Temp;
use Test::More tests => 4;

my $pkg = 'Genome::Model::Tools::Sx::Split::ByNs';
use_ok($pkg);

my $dir = File::Spec->join( Genome::Config::get('test_inputs'), 'Genome-Model-Tools-Sx', 'SplitByNs');

my $tmp_dir = Genome::Sys->create_temp_directory;

subtest 'fasta' => sub{
    plan tests => 4;

    my $in_fasta = File::Spec->join($dir, 'in.fasta');
    ok(-s $in_fasta, 'in fasta');

    my $out_fasta = File::Spec->join($tmp_dir, 'out.fasta');

    my $splitter = $pkg->execute(
        input  => [ $in_fasta ],
        output => [ $out_fasta ],
        number_of_ns => 10,
    );
    ok($splitter->result, 'execute');

    my $expected_fasta = File::Spec->join($dir, 'expected.fasta');
    ok(-s $expected_fasta, 'expected fasta');
    Genome::Utility::Test::compare_ok($out_fasta, $expected_fasta, 'output fasta matches expected');

};

subtest 'fastq' => sub{
    plan tests => 4;

    my $in_fastq = File::Spec->join($dir, 'in.fastq');
    ok(-s $in_fastq, 'in fastq');

    my $out_fastq = File::Spec->join($tmp_dir, 'out.fastq');

    my $splitter = $pkg->execute(
        input  => [ $in_fastq ],
        output => [ $out_fastq ],
        number_of_ns => 10,
    );
    ok($splitter->result, 'execute');

    my $expected_fastq = File::Spec->join($dir, 'expected.fastq');
    ok(-s $expected_fastq, 'expected fastq exists');
    Genome::Utility::Test::compare_ok($out_fastq, $expected_fastq, 'output fastq matches expected');

};

subtest 'buncha ns' => sub{
    plan tests => 4;

    my $in_fasta = File::Spec->join($dir, 'in.buncha-ns.fasta');
    ok(-s $in_fasta, 'in fasta');

    my $out_fasta = File::Spec->join($tmp_dir, 'out.bunch-ns.fasta');

    my $splitter = $pkg->execute(
        input  => [ $in_fasta ],
        output => [ $out_fasta ],
        number_of_ns => 50000, # 32767 is max for regex
    );
    ok($splitter->result, 'execute');

    my $expected_fasta = File::Spec->join($dir, 'expected.buncha-ns.fasta');
    ok(-s $expected_fasta, 'expected fasta exists');
    Genome::Utility::Test::compare_ok($out_fasta, $expected_fasta, 'output fasta matches expected');

};

done_testing();
