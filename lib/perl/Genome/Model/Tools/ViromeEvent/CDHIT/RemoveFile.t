#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok( 'Genome::Model::Tools::ViromeEvent::CDHIT::RemoveFile' );

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable';
ok( -d $data_dir, "Test suite data dir exists" );

my $rf = Genome::Model::Tools::ViromeEvent::CDHIT::RemoveFile->create(
    dir => $data_dir,
    );

ok( $rf, "Created virome screeing cd-hit remove files event" );

done_testing();

exit;
