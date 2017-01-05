#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Test::More tests => 5;

my $class = 'Genome::Config::AnalysisProject::Command::TakeOwnership';
use_ok($class);

my $ap = Genome::Config::AnalysisProject->create(
    name => 'Test Project',
    status => 'In Progress',
);
$ap->created_by('nobody');

my $cmd = $class->execute(analysis_projects => [$ap]);
ok($cmd->result, 'execute cmd');
is($ap->created_by, Genome::Sys->username, "AnP owner updated to current user");

my $cle_ap = Genome::Config::AnalysisProject->create(
    name => 'Test CLE Project',
    status => 'Pending',
    is_cle => 1,
);
$cle_ap->created_by('nobody');

my $cmd2 = $class->create(analysis_projects => [$cle_ap]);
isa_ok($cmd2, $class, 'created command');
my @errors = $cmd2->__errors__;
is(scalar(@errors), 1, 'found an error');

done_testing();
