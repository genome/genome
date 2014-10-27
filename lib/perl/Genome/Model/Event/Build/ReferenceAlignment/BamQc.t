#!/usr/bin/env genome-perl

use strict;
use warnings;

use version 0.77;
use above 'Genome';
use Test::More tests => 7;

use_ok('Genome::Model::Event::Build::ReferenceAlignment::BamQc');

is(Genome::Model::Event::Build::ReferenceAlignment::BamQc->_select_picard_version(1.13), Genome::Model::Tools::Picard->default_picard_version, "Picard version < 1.40 selected correctly.");
is(Genome::Model::Event::Build::ReferenceAlignment::BamQc->_select_picard_version(1.113), 1.113, "Picard version > 1.100 selected correctly.");
is(Genome::Model::Event::Build::ReferenceAlignment::BamQc->_bwa_mem_version_object("0.7.5a"), version->parse("v0.7.5_1"), "BWA mem lettered version munged correctly");
is(Genome::Model::Event::Build::ReferenceAlignment::BamQc->_bwa_mem_version_object("0.7.10"), version->parse("v0.7.10"), "BWA mem regular version munged correctly");
ok(Genome::Model::Event::Build::ReferenceAlignment::BamQc->_bwa_mem_version_object("0.7.10") > Genome::Model::Event::Build::ReferenceAlignment::BamQc->_bwa_mem_version_object("0.7.5a"), "0.7.10 is newer than 0.7.5a");
ok(Genome::Model::Event::Build::ReferenceAlignment::BamQc->_bwa_mem_version_object("0.7.5a") > Genome::Model::Event::Build::ReferenceAlignment::BamQc->_bwa_mem_version_object("0.7.5"), "0.7.5a is newer than 0.7.5");

