#!/usr/bin/env genome-perl

use strict;
use warnings;

use version 0.77;
use above 'Genome';
use Test::More tests => 5;

# This test was auto-generated because './Model/Event/Build/ReferenceAlignment/BamQc.pm'
# had no '.t' file beside it.  Please remove this test if you believe it was
# created unnecessarily.  This is a bare minimum test that just compiles Perl
# and the UR class.
use_ok('Genome::Model::Event::Build::ReferenceAlignment::BamQc');

is(Genome::Model::Event::Build::ReferenceAlignment::BamQc->_select_picard_version(1.13), Genome::Model::Tools::Picard->default_picard_version, "Picard version < 1.40 selected correctly.");
is(Genome::Model::Event::Build::ReferenceAlignment::BamQc->_select_picard_version(1.113), 1.113, "Picard version > 1.100 selected correctly.");
is(Genome::Model::Event::Build::ReferenceAlignment::BamQc->_bwa_mem_version_object("0.7.5a"), version->parse("v0.7.5_1"), "BWA mem lettered version munged correctly");
is(Genome::Model::Event::Build::ReferenceAlignment::BamQc->_bwa_mem_version_object("0.7.10"), version->parse("v0.7.10"), "BWA mem regular version munged correctly");

