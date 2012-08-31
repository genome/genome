#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

my $class = 'Genome::Model::Tools::TigraSv';

use_ok($class);

for my $ver ($class->available_tigrasv_versions) {
    ok(-x $class->path_for_tigrasv_version($ver), "tigra-sv $ver executable exists and is executable");
}

done_testing();
