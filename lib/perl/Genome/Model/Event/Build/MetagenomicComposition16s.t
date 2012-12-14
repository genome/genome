#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Event::Build::MetagenomicComposition16s') or die;
is(Genome::Model::Event::Build::MetagenomicComposition16s->bsub_rusage, "-R 'span[hosts=1]'", 'Busb rusage');

done_testing();
