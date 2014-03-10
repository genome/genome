#!/usr/bin/env perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 4;
use Genome::Command::Tester qw(run_and_diff);

run_and_diff(
    command => '$script_dir/converge/convergeStats.pl --model_group_id=66909 --outfile=$output_dir/out.tsv --verbose=1',
    results_version => '2014-01-30',
    eventual_class => 'Genome::Model::ClinSeq::Command::Converge::Stats',
);

