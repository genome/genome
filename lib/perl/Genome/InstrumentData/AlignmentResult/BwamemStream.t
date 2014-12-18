#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

# The actual test for this module is in Bwamem.t

my $pkg = 'Genome::InstrumentData::AlignmentResult::BwamemStream';
use_ok($pkg);


done_testing();
