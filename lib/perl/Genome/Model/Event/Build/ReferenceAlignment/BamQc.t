#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 1;

# This test was auto-generated because './Model/Event/Build/ReferenceAlignment/BamQc.pm'
# had no '.t' file beside it.  Please remove this test if you believe it was
# created unnecessarily.  This is a bare minimum test that just compiles Perl
# and the UR class.
use_ok('Genome::Model::Event::Build::ReferenceAlignment::BamQc');

is(Genome::Model::Event::Build::ReferenceAlignment::BamQc->_select_picard_version(1.13), Genome::Model::Tools::Picard->default_picard_version, "Picard version < 1.40 selected correctly.");
is(Genome::Model::Event::Build::ReferenceAlignment::BamQc->_select_picard_version(1.113), 1.113, "Picard version > 1.100 selected correctly.");
