#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Test::More tests => 6;

my $class = 'Genome::Config::AnalysisProject::Command::Uncomplete';
use_ok($class);

my $ap = Genome::Config::AnalysisProject->create(
    name => 'Test Project',
    status => 'Completed',
);
is($ap->status, 'Completed', "AnP status is 'In Progress'");

my $cmd = $class->execute(analysis_projects => [$ap]);
ok($cmd->result, 'execute cmd');
my $expected_status = 'In Progress';
is($ap->status, $expected_status, "AnP status set to '$expected_status'");

my $other_ap = Genome::Config::AnalysisProject->create(
    name => 'Test Project',
    status => 'Completed',
);
$other_ap->created_by('nobody');

my $cmd2 = $class->create(analysis_projects => [$other_ap]);
isa_ok($cmd2, $class, 'created command');
my @errors = $cmd2->__errors__;
is(scalar(@errors), 1, 'found an error');

done_testing();
