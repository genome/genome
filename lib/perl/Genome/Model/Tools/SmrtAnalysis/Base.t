#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

unless (Genome::Sys->username eq 'smrtanalysis') {
  plan skip_all => "this test is only runnable by user smrtanalysis"
}
plan tests => 1;

use_ok('Genome::Model::Tools::SmrtAnalysis::Base');
