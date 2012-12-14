#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Compare;
require File::Temp;
use Test::More;

use_ok('Genome::Model::Tools::Sx::Metrics::Basic') or die;

# Dir/files
my $testdir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/Metrics/Basic/v1';
my $example_file = $testdir.'/metrics.txt';
my $legacy_file = $testdir.'/metrics.legacy';

# Create, add sequences
my $metrics = Genome::Model::Tools::Sx::Metrics::Basic->create(id => 'TEST');
ok($metrics, 'create');
ok($metrics->add_sequence({seq => 'AAGGCCTT',}), 'add seq');
is($metrics->count, 1, 'count');
is($metrics->bases, 8, 'bases');
ok($metrics->add_sequences([{seq => 'AAGGCCTT',}]), 'add seqs');
is($metrics->count, 2, 'count');
is($metrics->bases, 16, 'bases');
# Save to file
my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $file = $tmpdir.'/metrics';
ok($metrics->to_file($file), 'to file');
ok(-s $file, 'metrics file exists');
is(File::Compare::compare($file, $example_file), 0, 'metrics file matches');

# Reload the metrics from the file
my $metrics_from_file = Genome::Model::Tools::Sx::Metrics->from_file($file);
ok($metrics_from_file, 'get metrics from file');
is($metrics_from_file->class, $metrics->class, 'metrics from file class is correct');
is_deeply($metrics_from_file, $metrics, 'metrics from file match');
is_deeply(
    { id => $metrics_from_file->id, bases => $metrics_from_file->bases, count => $metrics_from_file->count },
    { id => 'TEST', bases => 16, count => 2 }, 
    'metrics from file match',
);

# Legacy from file
my $metrics_from_file_legacy = Genome::Model::Tools::Sx::Metrics->from_file($legacy_file);
ok($metrics_from_file_legacy, 'get metrics from legacy file');
is($metrics_from_file->class, $metrics->class, 'metrics from legacy file class is correct');
isnt($metrics_from_file_legacy->id, 'TEST', 'Id is not TEST');
is_deeply({ bases => $metrics->bases, count => $metrics->count }, { bases => 16, count => 2 }, 'metrics from legacy file match');

#print "$tmpdir\n"; <STDIN>;
done_testing();
