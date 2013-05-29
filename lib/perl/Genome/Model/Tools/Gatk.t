#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 4;

# This test was auto-generated because './Model/Tools/Gatk.pm'
# had no '.t' file beside it.  Please remove this test if you believe it was
# created unnecessarily.  This is a bare minimum test that just compiles Perl
# and the UR class.
use_ok('Genome::Model::Tools::Gatk');

ok(Genome::Model::Tools::Gatk->is_legacy_version("v1"), "v1 is a legacy version");
ok(Genome::Model::Tools::Gatk->is_legacy_version("5777"), "5777 is a legacy version");
ok(!Genome::Model::Tools::Gatk->is_legacy_version("2.4"), "2.4 is not a legacy version");
