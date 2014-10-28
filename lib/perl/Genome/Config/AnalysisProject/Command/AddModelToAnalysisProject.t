#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Test::More tests => 29;

use Genome::Test::Factory::AnalysisProject;
use Genome::Test::Factory::Model::ReferenceAlignment;

my $class = 'Genome::Config::AnalysisProject::Command::AddModelToAnalysisProject';

use_ok($class);

my $analysis_project = Genome::Test::Factory::AnalysisProject->setup_object(status => 'Pending');
my ($profile_item) = $analysis_project->config_items;
$profile_item->status('disabled');
my @models = map Genome::Test::Factory::Model::ReferenceAlignment->setup_object, (1..3);

my $cmd1 = $class->create(
    profile_item => $profile_item,
    models => \@models,
);
isa_ok($cmd1, $class, 'created command');
ok(!$cmd1->execute, 'command fails on pending analysis project');
for my $m (@models) {
    ok(!$m->analysis_projects, 'no analysis project assigned to models');
}

my $cmd2 = $class->create(
    profile_item => $profile_item,
    models => \@models,
    allow_projects_not_in_progress => 1,
);
isa_ok($cmd2, $class, 'created command');
ok($cmd2->execute, 'command succeeds when pending analysis project allowed');
for my $m (@models) {
    is($m->analysis_projects, $analysis_project, 'analysis project assigned');
}

my $cmd3 = $class->create(
    profile_item => $profile_item,
    models => \@models,
    allow_projects_not_in_progress => 1,
);
isa_ok($cmd3, $class, 'created command');
ok($cmd3->execute, 'command succeeds when models already assigned to that project');


my $other_analysis_project = Genome::Test::Factory::AnalysisProject->setup_object(status => 'In Progress');
my ($other_profile_item) = $other_analysis_project->config_items;
$other_profile_item->status('disabled');
my $other_model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object();

my $cmd4 = $class->create(
    profile_item => $other_profile_item,
    models => [@models, $other_model],
);
isa_ok($cmd4, $class, 'created command');
ok(!$cmd4->execute, 'command fails when models already assigned to a different project');
ok(!$other_model->analysis_projects, 'unassigned model not assigned to');
for my $m (@models) {
    is($m->analysis_projects, $analysis_project, 'analysis project remains properly assigned');
    is($m->config_profile_items, $profile_item, 'model linked to correct configuration');
}

$other_profile_item->status('active');
my $cmd5 = $class->create(
    profile_item => $other_profile_item,
    models => [$other_model],
);
isa_ok($cmd5, $class, 'created command');
ok(!$cmd5->execute, 'command fails when attempting to assign to an active config profile item');
ok(!$other_model->analysis_projects, 'unassigned model not assigned to');

$other_profile_item->status('disabled');
my $cmd6 = $class->create(
    profile_item => $other_profile_item,
    models => [$other_model],
);
isa_ok($cmd6, $class, 'created command');
ok($cmd6->execute, 'command succeeds in "normal" case');
is($other_model->analysis_projects, $other_analysis_project, 'model assigned appropriately');
is($other_model->config_profile_items, $other_profile_item, 'model linked to correct configuration');
