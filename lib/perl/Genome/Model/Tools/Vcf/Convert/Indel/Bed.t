#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

use above 'Genome';

use_ok('Genome::Model::Tools::Vcf::Convert::Indel::Bed');

done_testing();
