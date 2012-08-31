#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
    
use Test::More;

use_ok('Genome::Report::FromLegacy');

my $generator = Genome::Report::FromLegacy->create(
    properties_file => $ENV{GENOME_TEST_INPUTS} . '/Genome-Report-FromLegacy/Legacy_Report/properties.stor',
);
ok($generator, 'create legacy generator');

my $report = $generator->generate_report;
ok($report, 'generated report');

done_testing();
exit;

