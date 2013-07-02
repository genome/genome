#!/usr/bin/env perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 4;
use Genome::Command::Tester qw(run_and_diff);

run_and_diff(
    #command => 'Genome::Model::SomaticVariation::Command::Compare::Variants',
    results_version => '2013-07-02',
);

