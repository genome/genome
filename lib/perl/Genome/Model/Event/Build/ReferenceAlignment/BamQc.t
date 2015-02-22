#!/usr/bin/env genome-perl

use strict;
use warnings;

use version 0.77;
use above 'Genome';
use Genome::Test::Factory::ProcessingProfile::ReferenceAlignment;
use Test::More;

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

my $test_pp = Genome::Test::Factory::ProcessingProfile::ReferenceAlignment->setup_object(
);

for my $aligner (qw(bwamem bwamem-stream)) {
    $test_pp->read_aligner_name($aligner);
    $test_pp->read_aligner_version('0.7.5');

    is(
        $bamqc_class->_select_error_rate_version_for_pp($test_pp),
        Genome::Model::Tools::BioSamtools::ErrorRate->default_errorrate_version(),
        "Default ErrorRate version chosen correctly for aligner $aligner"
    );

    $test_pp->read_aligner_version('0.7.6');


    is($bamqc_class->_should_skip_bam_qc($test_pp),
        0,
        'BamQc correctly chosen to run'
    );

    is($bamqc_class->_select_error_rate_version_for_pp($test_pp),
        '1.0a2',
        'New ErrorRate version chosen correctly'
    );
}

$test_pp->read_aligner_version();
$test_pp->read_aligner_name('bsmap');

is($bamqc_class->_should_skip_bam_qc($test_pp),
    1,
    'BamQc correctly skipped for bsmap'
);

$test_pp->read_aligner_name('imported');

is($bamqc_class->_should_skip_bam_qc($test_pp),
    1,
    'BamQc correctly skipped for imported'
);

done_testing();
