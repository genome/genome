#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_NO_REQUIRE_USER_VERIFY} = 1;
}

use strict;
use warnings;

use above "Genome";
use Test::More;

use Genome::Test::Factory::AnalysisProject;

use_ok("Genome::Model::Build::Command::AbandonAndQueue") or die;

class Genome::Model::Tester { is => 'Genome::ModelDeprecated', };

my $s = Genome::Sample->create(name => 'TEST-' . __FILE__ . "-$$");
ok($s, "made a test sample");

my $p = Genome::ProcessingProfile::Tester->create(
    name => 'Tester Test for Testing',
);
ok($p, "made a test processing profile");

my $m = Genome::Model::Tester->create(
    processing_profile_id => $p->id,
    subject_class_name => ref($s),
    subject_id => $s->id,
    build_requested => 0,
);
ok($m, "made a test model");
ok(!$m->build_requested, 'build is not requested');

my $anp = Genome::Test::Factory::AnalysisProject->setup_object();
$anp->add_model_bridge(model_id => $m->id);

my $b1 = $m->add_build();
ok($b1, "made test build 1");

my $exit_code1 = eval { 
    Genome::Model::Build::Command::AbandonAndQueue->_execute_with_shell_params_and_return_exit_code('--reason', 'to test the abandon-and-queue command', '--', $b1->id);
};
ok(!$@, "the command did not crash");
is($exit_code1, 0, "command believes it succeeded");
ok($m->build_requested, 'build is requested');

done_testing();
