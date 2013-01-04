#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use above "Genome";
use Genome::Model::Tools::Annotate;

# Using Genome::Model::Tools::Annotate should not enable cache pruning options
# because the class is loaded as part of calculating sub-commands.
is(UR::Context->object_cache_size_lowwater, undef, 'lowwater cache size is not set');
is(UR::Context->object_cache_size_highwater, undef, 'highwater cache size is not set');

my ($low, $high) = Genome::Model::Tools::Annotate->object_cache_sizes();
is($low, undef, 'lowwater cache size is not set');
is($high, undef, 'highwater cache size is not set');

Genome::Model::Tools::Annotate->object_cache_sizes(10, 100);
is(UR::Context->object_cache_size_lowwater, 10, 'lowwater cache size set to 10');
is(UR::Context->object_cache_size_highwater, 100, 'highwater cache size set to 100');

Genome::Model::Tools::Annotate->object_cache_sizes(undef, undef);
is(UR::Context->object_cache_size_lowwater, undef, 'lowwater cache size unset');
is(UR::Context->object_cache_size_highwater, undef, 'highwater cache size unset');

done_testing();
