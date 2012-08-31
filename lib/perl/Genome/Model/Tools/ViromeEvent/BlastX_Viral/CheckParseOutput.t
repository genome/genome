#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::ViromeEvent::BlastX_Viral::CheckParseOutput');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17/Titanium17_undecodable';
ok( -d $data_dir, "Test suite data dir exists" );

ok( -s '/gscmnt/sata835/info/medseq/virome/taxonomy_db', "Taxonomy db exists" );

#create
my $c = Genome::Model::Tools::ViromeEvent::BlastX_Viral::CheckParseOutput->create(
    dir => $data_dir,
    );
ok($c, "Created blastx viral check parse output event");

done_testing();

exit;
