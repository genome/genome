#!/usr/bin/env genome-perl

use strict;
use warnings;

use version 0.77;
use above 'Genome';
use Test::More tests => 7;

my $bamqc_class = 'Genome::Model::Event::Build::ReferenceAlignment::BamQc';

use_ok($bamqc_class);

is(
    $bamqc_class->_select_picard_version(1.13),
    Genome::Model::Tools::Picard->default_picard_version,
    'Picard version < 1.40 selected correctly.'
);

is(
    $bamqc_class->_select_picard_version(1.113),
    1.113,
    'Picard version > 1.100 selected correctly.'
);

is(
    $bamqc_class->_bwa_mem_version_object('0.7.5a'),
    version->parse('v0.7.5_1'),
    'BWA mem lettered version munged correctly'
);

is(
    $bamqc_class->_bwa_mem_version_object('0.7.10'),
    version->parse('v0.7.10'),
    'BWA mem regular version munged correctly'
);

ok(
    $bamqc_class->_bwa_mem_version_object('0.7.10') > $bamqc_class->_bwa_mem_version_object('0.7.5a'),
    '0.7.10 is newer than 0.7.5a'
);

ok(
    $bamqc_class->_bwa_mem_version_object('0.7.5a') > $bamqc_class->_bwa_mem_version_object('0.7.5'),
    '0.7.5a is newer than 0.7.5'
);

