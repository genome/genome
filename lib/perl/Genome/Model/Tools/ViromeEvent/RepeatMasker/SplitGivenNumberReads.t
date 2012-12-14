#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::ViromeEvent::RepeatMasker::SplitGivenNumberReads');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable';
ok( -d $data_dir, "Test suite data dir exists" );

my $c = Genome::Model::Tools::ViromeEvent::RepeatMasker::SplitGivenNumberReads->create(
    dir => $data_dir,
    );
ok($c, "Created repeat masker split given number of reads event");

done_testing();
