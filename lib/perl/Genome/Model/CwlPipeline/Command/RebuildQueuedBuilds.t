#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 1;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

my $class = 'Genome::Model::Command::Services::BuildQueuedModels';
use_ok($class);
