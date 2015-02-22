#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper 'Dumper';
use Test::More;

use_ok('Genome::WorkOrder');

my $wo = Genome::WorkOrder->get(2231383);
ok($wo, 'Got work order');

my @woi = $wo->items;
ok(@woi, 'Got work order items');

done_testing();
