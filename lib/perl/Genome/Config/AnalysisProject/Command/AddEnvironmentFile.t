#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Test::More tests => 3;
use Test::Exception;

my $class = 'Genome::Config::AnalysisProject::Command::AddEnvironmentFile';
use_ok($class);

my $ap = Genome::Config::AnalysisProject->create(
    name => 'Test Project',
    id => '-500',
);

my $file = Genome::Sys->create_temp_file_path();
Genome::Sys->write_file($file, 'no_commit: 1');

my @params = (
    analysis_project => $ap,
    environment_file => $file,
);

my $cmd = $class->create(@params);

$cmd->execute();

ok($ap->environment_config_dir, 'environment file was put in place');

my $cmd2 = $class->create(@params);

throws_ok(sub {$cmd2->execute}, qr/already has an existing environment file/, 'cannot add environment file a second time');

done_testing();

1;
