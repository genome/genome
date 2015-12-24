#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Test::More tests => 5;

my $class = 'Genome::Config::AnalysisProject::Command::Complete';
use_ok($class);

my $ap = Genome::Config::AnalysisProject->create(
    name => 'Test Project',
    status => 'In Progress',
);
is($ap->status, 'In Progress', "AnP status is 'In Progress'");

my $cmd = $class->execute(analysis_projects => [$ap]);
ok($cmd->result, 'execute cmd');
my $expected_status = 'Completed';
is($ap->status, $expected_status, "AnP status set to '$expected_status'");
is($cmd->status_message, 'Completed: '.$ap->__display_name__, 'Status message about releasing');

done_testing();
