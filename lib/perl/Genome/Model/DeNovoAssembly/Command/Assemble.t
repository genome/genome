#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

diag('TESTS ARE DONE IN BUILD SUBCLASSES');
use_ok('Genome::Model::DeNovoAssembly::Command::Assemble') or die;

done_testing();
