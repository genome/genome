#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::ViromeEvent::BlastN::CheckParseOutput');

ok( -s '/gscmnt/sata835/info/medseq/virome/taxonomy_db', "Taxonomy blast db exists" );

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable/';
ok ( -d $data_dir, "Test suite data dir exists" );

my $c = Genome::Model::Tools::ViromeEvent::BlastN::CheckParseOutput->create(
    dir => $data_dir,
    );
ok($c, "Created virome event blast-n check-parse-output event");

done_testing();

exit;
