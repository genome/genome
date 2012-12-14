#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 2;

use above 'Genome';

use_ok('Genome::Capture::Region');
my $region = Genome::Capture::Region->get(77373);
isa_ok($region,'Genome::Capture::Region');
