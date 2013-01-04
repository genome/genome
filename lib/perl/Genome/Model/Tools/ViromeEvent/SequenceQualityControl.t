#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::ViromeEvent::SequenceQualityControl');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17';
ok( -d $data_dir, "Test suite data dir exists" );

my $c = Genome::Model::Tools::ViromeEvent::SequenceQualityControl->create(
    dir => '/gscmnt/sata835/info/medseq/virome/test17',
    );

ok( $c, "Created sequence-quality-control event" );

done_testing();
