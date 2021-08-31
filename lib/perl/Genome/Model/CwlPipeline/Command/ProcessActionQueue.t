#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 6;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';
use Genome::Test::Factory::AnalysisProject;
use Genome::Test::Factory::Build;
use Genome::Test::Factory::Model::CwlPipeline;


my $class = 'Genome::Model::CwlPipeline::Command::ProcessActionQueue';
use_ok($class);

my $anp = Genome::Test::Factory::AnalysisProject->setup_object();

my $model = Genome::Test::Factory::Model::CwlPipeline->setup_object();
my $bridge = Genome::Config::AnalysisProject::ModelBridge->create(
    model_id => $model->id,
    analysis_project_id => $anp->id,
);

my $build = Genome::Test::Factory::Build->setup_object( model_id => $model->id, status => 'Failed');

$build->action_requested('abandon');
is($build->action_requested, 'abandon', 'set up build to abandon');

my $cmd = $class->create(
    submit_jobs => 0,
    max_count => 5,
);
isa_ok($cmd, $class, 'created command');

#avoid finding production objects
my $original_quc = UR::Context->query_underlying_context();
UR::Context->query_underlying_context(0);

ok($cmd->execute, 'executed command');

is($build->action_requested, 0, 'action requested cleared');
is($build->status, 'Abandoned', 'build was abandoned');

done_testing();
