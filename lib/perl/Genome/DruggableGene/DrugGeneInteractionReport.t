#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

use above "Genome";

my $class = 'Genome::DruggableGene::DrugGeneInteractionReport';
use_ok($class);

done_testing();
