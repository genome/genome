#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Test::More tests => 1;

my $class = 'Genome::Config::AnalysisMenu::Command::Item::Create';
use_ok($class);

1;
