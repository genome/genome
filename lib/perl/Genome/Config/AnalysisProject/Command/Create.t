#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
    $ENV{UR_NO_REQUIRE_USER_VERIFY} = 1;
};

use Test::More;

use above 'Genome';

use_ok('Genome::Config::AnalysisProject::Command::Create');

my $cmd = Genome::Config::AnalysisProject::Command::Create->create(
    name => 'test proj',
    is_production => 1,
);
ok($cmd, 'constructed create command');
isa_ok($cmd, 'Genome::Config::AnalysisProject::Command::Create');
my $res = $cmd->execute;
ok($res, 'command executed successfully');
isa_ok($res, 'Genome::Config::AnalysisProject', 'command returned a Genome::Config::AnalysisProject');

my $cmd2 = Genome::Config::AnalysisProject::Command::Create->create(
    name => 'test proj non-production',
    is_production => 0,
);
ok($cmd2, 'constructed second create command');
isa_ok($cmd2, 'Genome::Config::AnalysisProject::Command::Create');
my $res2 = $cmd2->execute;
ok($res2, 'second command executed successfully');
isa_ok($res2, 'Genome::Config::AnalysisProject', 'command returned a Genome::Config::AnalysisProject');
is($res2->run_as, Genome::Sys->username, 'run_as set correctly');

my $cmd3 = Genome::Config::AnalysisProject::Command::Create->create(
    name => 'test bad create options',
    is_production => 1,
    is_cle => 1,
);
ok($cmd3, 'constructed third create command');
isa_ok($cmd3, 'Genome::Config::AnalysisProject::Command::Create');
my $res3 = eval { $cmd3->execute };
my $err3 = $@;
ok(!$res3, 'third command fails with bad options');
ok($err3, 'error thrown with bad options');

done_testing();
