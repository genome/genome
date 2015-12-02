#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More tests => 1;

diag('Tests are in the inputs factory!');
use_ok('Genome::InstrumentData::Command::Import::Inputs::Factory');

done_testing();
