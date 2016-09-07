#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;   
};

use strict;
use warnings;

use above 'Genome';
use Test::More;

use Genome::Test::Factory::AnalysisProject;

use_ok('Genome::Model::Build::Command::Queue') or die;

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
);
ok($m, "made a test model");
my $model_id = $m->id;

my $anp = Genome::Test::Factory::AnalysisProject->setup_object();
$anp->add_model_bridge(model_id => $model_id);

my $reason = 'testing the queue command';
my $cmd = Genome::Model::Build::Command::Queue->execute(
    models => [$m],
    reason => $reason,
);
ok($cmd->result, "command believes it succeeded");

is($m->build_requested, 1, 'The command requested a build for the model');
my @n = $m->notes;
my ($note) = grep($_->header_text eq 'build_requested', @n);
ok($note, 'found a note about the build being requested');
is($note->body_text, $reason, 'note has expected reason');

done_testing();
