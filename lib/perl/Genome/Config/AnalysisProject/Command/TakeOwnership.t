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
    status => 'Hold',
    created_by => 'foo',
);
is($ap->created_by, 'foo', "AnP created_by is 'foo'");

my $cmd = $class->execute(analysis_projects => [$ap]);
ok($cmd->result, 'execute cmd');
my $expected_user = $ENV{USER};
is($ap->created_by, $expected_user, "AnP created_by set to '$expected_user'");
is($cmd->status_message, 'Updated created_by from foo to '. $expected_user .' for project '.$ap->__display_name__);     

done_testing();
