#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Report::Table') or die;

my $generator = Genome::Model::Report::Table->create(
    name => 'Table Test',
    description => 'A testing of the table report generator',
    headers => [qw/ model_id build_id status date /],
    row_name => 'build',
    rows => [
        [qw| 2816929867 98421139 Succeeded 2009-08-27 |],
        [qw| 2816929867 98421140 Succeeded 2009-08-28 |],
        [qw| 2816929867 98421141 Succeeded 2009-12-29 |],
        [qw| 2816929868 98421142 Succeeded 2009-08-27 |],
        [qw| 2816929868 98421143 Abandoned 2009-08-30 |],
    ],
);
ok($generator, 'create');
my $report = $generator->generate_report;
ok($report, 'generate report');

done_testing();
