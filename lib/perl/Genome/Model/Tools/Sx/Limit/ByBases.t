#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Compare;
use Test::More;

# use
use_ok('Genome::Model::Tools::Sx::Limit::ByBases') or die;

# fail
ok(!Genome::Model::Tools::Sx::Limit::ByBases->execute(), 'failed w/o bases');
ok(!Genome::Model::Tools::Sx::Limit::ByBases->execute(bases => 'all'), 'failed  w/ bases => all');
ok(!Genome::Model::Tools::Sx::Limit::ByBases->execute(bases => 0), 'failed w/ bases => 0');

# Files
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/LimitByBases/v2';
my $in_fastq = $dir.'/in.fastq';
ok(-s $in_fastq, 'in fastq');
my $example_fastq = $dir.'/out.fastq';
ok(-s $example_fastq, 'example fastq');
my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);

# Ok - limit regular
my $out_fastq = $tmp_dir.'/out.fastq';
my $out_metrics = $tmp_dir.'/out.metrics';
my $limiter = Genome::Model::Tools::Sx::Limit::ByBases->create(
    input  => [ $in_fastq ],
    output => [ $out_fastq ],
    output_metrics => $out_metrics,
    bases => 6000, # 80% in.fastq has 100 sequences w/ 75 bases each for 7500 total
);
ok($limiter, 'create limiter');
isa_ok($limiter, 'Genome::Model::Tools::Sx::Limit::ByBases');
ok($limiter->execute, 'execute limiter');
is(File::Compare::compare($example_fastq, $out_fastq), 0, "fastq limited as expected");
my $metrics = Genome::Model::Tools::Sx::Metrics->from_file($out_metrics);
is_deeply([$metrics->bases, $metrics->count], [qw/ 6000 80 /], 'randomly limited the correct number of sequences and bases');

# Ok - random
my $out_random_fastq = $tmp_dir.'/out.random.fastq';
my $out_random_metrics = $tmp_dir.'/out.random.metrics';
$limiter = Genome::Model::Tools::Sx::Limit::ByBases->create(
    input  => [ $in_fastq ],
    output => [ $out_random_fastq ],
    output_metrics => $out_random_metrics,
    bases => 375, # 5% in.fastq has 100 sequences w/ 75 bases each for 7500 total
    incoming_sequences => $dir.'/in.metrics',
);
ok($limiter, 'create limiter');
isa_ok($limiter, 'Genome::Model::Tools::Sx::Limit::ByBases');
ok($limiter->execute, 'execute limiter');
ok(-s $out_random_fastq, 'created random fastq');
is(File::Compare::compare($out_random_fastq, $out_fastq), 1, "random fastq is different than example");
$metrics = Genome::Model::Tools::Sx::Metrics->from_file($out_random_metrics);
cmp_ok($metrics->bases, '<=', 375, 'randomly limited the correct number of bases');
cmp_ok($metrics->count, '<=', 5, 'randomly limited the correct number of sequences');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
