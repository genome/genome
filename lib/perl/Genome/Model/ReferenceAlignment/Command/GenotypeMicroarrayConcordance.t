#!/usr/bin/env genome-perl

use strict;
use warnings 'FATAL';

use above 'Genome';

use Test::More;
plan tests => 1;

my $class = 'Genome::Model::ReferenceAlignment::Command::GenotypeMicroarrayConcordance';
use_ok($class) or die;

done_testing();
