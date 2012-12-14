#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 6;

use above 'Genome';

use_ok('Genome::Capture::Target');
my $target = Genome::Capture::Target->get(403382);
isa_ok($target,'Genome::Capture::Target');
is($target->region_id,'77373','Found region id');
isa_ok($target->region,'Genome::Capture::Region');
is($target->tag_id,'2764758292','Found tag id ');
is($target->pse_id,'95372894','Found pse id');
