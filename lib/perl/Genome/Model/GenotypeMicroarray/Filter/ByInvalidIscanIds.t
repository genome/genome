#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::GenotypeMicroarray::Filter::ByInvalidIscanIds') or die;

my $filter = Genome::Model::GenotypeMicroarray::Filter::ByInvalidIscanIds->create();
ok($filter, 'create filter');
ok($filter->filter({id => 'nathan'}), 'did not filter id not in list');
ok(!$filter->filter({id => 'rs1010408'}), 'filter id in list');

done_testing();
