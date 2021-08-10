#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
}

use above 'Genome';

use Test::More tests => 1;

use_ok('Genome::Disk::Group') or die;

done_testing();
