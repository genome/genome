#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Test::More tests => 5;

use File::Spec;
use Genome::Test::Factory::DiskAllocation;

my $class = 'Genome::Config::AnalysisProject::Command::Release';
use_ok($class);

my $ap = Genome::Config::AnalysisProject->create(
    name => 'Test Project',
    status => 'Hold',
);

use Genome::Config;
my $da = Genome::Test::Factory::DiskAllocation->setup_object(owner => $ap);
Genome::Sys->create_directory(File::Spec->join($da->absolute_path, Genome::Config::config_subpath));

is($ap->status, 'Hold', "AnP status is 'Hold'");

my $cmd = $class->execute(analysis_projects => [$ap]);
ok($cmd->result, 'execute cmd');
my $expected_status = 'In Progress';
is($ap->status, $expected_status, "AnP status set to '$expected_status'");
is($cmd->status_message, 'Released: '.$ap->__display_name__, 'Status message about releasing');

done_testing();
