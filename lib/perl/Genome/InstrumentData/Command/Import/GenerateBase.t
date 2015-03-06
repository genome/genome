#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';
use Test::More;

# Tests are in subclasses
use_ok('Genome::InstrumentData::Command::Import::GenerateBase') or die;

done_testing();
