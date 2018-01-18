#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';

use Genome::Model::CwlPipeline::Command::Run; #to override later
use Genome::Test::Factory::AnalysisProject;
use Genome::Test::Factory::Model::CwlPipeline;

use Sub::Install;
use Test::More tests => 4;

my $class = 'Genome::Model::CwlPipeline::Command::Restart';
use_ok($class);

my $anp = Genome::Test::Factory::AnalysisProject->setup_object;
my $model = Genome::Test::Factory::Model::CwlPipeline->setup_object(
    name => 'testing restart command',
);
$model->add_analysis_project_bridge(analysis_project_id => $anp->id);
my $build = Genome::Model::Build::CwlPipeline->create(
    model_id => $model->id,
);

Genome::Report::Email->silent;
my $guard = Genome::Config::set_env('workflow_builder_backend', 'inline');
my $override = Sub::Install::reinstall_sub({
    into => 'Genome::Model::CwlPipeline::Command::Run',
    as => '_execute_body',
    code => sub { return; }, #let's fail.
});

$build->start;
$build->status('Failed'); #inline test ends up unstartable, but this is atypical.

undef $override;
$override = Sub::Install::reinstall_sub({
    into => 'Genome::Model::CwlPipeline::Command::Run',
    as => '_execute_body',
    code => sub { return 1; }, #let's succeed.
});
my $cmd = $class->create(
    builds => [$build],
);
isa_ok($cmd, $class, 'created command');
ok($cmd->execute, 'executed command, and it succeeded');

is($build->status, 'Succeeded', 'build succeeded on second attempt');
