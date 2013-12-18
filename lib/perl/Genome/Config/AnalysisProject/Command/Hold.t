#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Test::More;

my $class = 'Genome::Config::AnalysisProject::Command::Hold';
use_ok($class);

my $ap = Genome::Config::AnalysisProject->create(
    name => 'Test Project'
);

is($ap->status, 'Pending', 'initial status should be pending');

my $cmd = $class->create(analysis_projects => [$ap]);
$cmd->execute();

is($ap->status, 'Hold', 'it should set the status to hold');

done_testing();

1;
