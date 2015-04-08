#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

my $class = 'Genome::Model::DeNovoAssembly::Build::ProcessInstrumentDataBase';
use_ok($class) or die;

my $lsf_resource = $class->lsf_resource;
ok($lsf_resource, 'lsf resource');

done_testing();
