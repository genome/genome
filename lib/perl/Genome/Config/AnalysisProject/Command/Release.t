#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Test::More;

my $class = 'Genome::Config::AnalysisProject::Command::Release';
use_ok($class);

my $ap = Genome::Config::AnalysisProject->create(
    name => 'Test Project',
    status => 'Hold',
);

is($ap->status, 'Hold', 'initial status should be pending');

my $cmd = $class->create(analysis_projects => [$ap]);
$cmd->execute();

is($ap->status, 'In Progress', 'it should set the status to in progress');

done_testing();

1;
