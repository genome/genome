#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Compare;
use Test::More;

# use
use_ok('Genome::Model::Tools::Sx::Limit::ByCount') or die;

# fail
ok(!Genome::Model::Tools::Sx::Limit::ByCount->execute(), 'failed w/o count');
ok(!Genome::Model::Tools::Sx::Limit::ByCount->execute(count => 'all'), 'failed  w/ count => all');
ok(!Genome::Model::Tools::Sx::Limit::ByCount->execute(count => 0), 'failed w/ count => 0');

# Files
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx';
my $in_fastq = $dir.'/in.fastq';
ok(-s $in_fastq, 'in fastq');
my $example_fastq = $dir.'/limit_by_coverage.example.fastq';
ok(-s $example_fastq, 'example fastq');
my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);

# Ok - count
my $out_fastq = $tmp_dir.'/out.count.fastq';
my $limiter = Genome::Model::Tools::Sx::Limit::ByCount->create(
    input  => [ $in_fastq ],
    output => [ $out_fastq ],
    count => 1234,
);
ok($limiter, 'create limiter');
isa_ok($limiter, 'Genome::Model::Tools::Sx::Limit::ByCount');
ok($limiter->execute, 'execute limiter');
is(File::Compare::compare($example_fastq, $out_fastq), 0, "fastq limited as expected");

#print "$tmp_dir\n"; <STDIN>;
done_testing();
exit;

