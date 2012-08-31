#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::MetagenomicCompositionShotgun::Command::TaxonomyReport') or die;

# fails
ok(!Genome::Model::MetagenomicCompositionShotgun::Command::TaxonomyReport->execute(), 'failed to execute w/o params');

done_testing();
exit;

