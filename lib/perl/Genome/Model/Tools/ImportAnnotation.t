#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use above "Genome";
use Genome::Model::Tools::ImportAnnotation;

# Using Genome::Model::Tools::ImportAnnotation should not enable cache pruning options
# because the class is loaded as part of calculating sub-commands.
ok(!UR::Context->object_cache_size_highwater, 'object_cache_size_highwater is not set');

done_testing();

