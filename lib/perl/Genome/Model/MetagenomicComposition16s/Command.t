#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper;
use Genome::Model::Build::MetagenomicComposition16s::TestBuildFactory;
use Test::More;

use_ok('Genome::Model::MetagenomicComposition16s::Command');

done_testing();
