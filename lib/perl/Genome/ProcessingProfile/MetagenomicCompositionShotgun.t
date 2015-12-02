#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 1;

# This model type is deprecated
use_ok('Genome::ProcessingProfile::MetagenomicCompositionShotgun');
