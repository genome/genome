#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Test::More tests => 30;

use Genome::Test::Factory::AnalysisProject;
use Genome::Test::Factory::Model::ReferenceAlignment;

my $class = 'Genome::Config::AnalysisProject::Command::AddModel';

use_ok($class);

my $analysis_project = Genome::Test::Factory::AnalysisProject->setup_object(status => 'Pending');
my $profile_item = add_config($analysis_project, 'Genome::Model::ReferenceAlignment');
my @models = map Genome::Test::Factory::Model::ReferenceAlignment->setup_object, (1..3);

my $cmd2 = $class->create(
    profile_item => $profile_item,
    models => \@models,
);
isa_ok($cmd2, $class, 'created command');
ok($cmd2->execute, 'command succeeds when pending analysis project');
for my $m (@models) {
    is($m->analysis_project, $analysis_project, 'analysis project assigned');
}

my $cmd3 = $class->create(
    profile_item => $profile_item,
    models => \@models,
);
isa_ok($cmd3, $class, 'created command');
ok($cmd3->execute, 'command succeeds when models already assigned to that project');


my $other_analysis_project = Genome::Test::Factory::AnalysisProject->setup_object(status => 'In Progress');
my $other_profile_item = add_config($other_analysis_project, 'Genome::Model::ReferenceAlignment');
my $other_model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object();

my $cmd4 = $class->create(
    profile_item => $other_profile_item,
    models => [@models, $other_model],
);
isa_ok($cmd4, $class, 'created command');
ok(!$cmd4->execute, 'command fails when models already assigned to a different project');
ok(!$other_model->analysis_project, 'unassigned model not assigned to');
for my $m (@models) {
    is($m->analysis_project, $analysis_project, 'analysis project remains properly assigned');
    is($m->config_profile_items, $profile_item, 'model linked to correct configuration');
}

$other_profile_item->status('active');
my $cmd5 = $class->create(
    profile_item => $other_profile_item,
    models => [$other_model],
);
isa_ok($cmd5, $class, 'created command');
ok(!$cmd5->execute, 'command fails when attempting to assign to an active config profile item');
ok(!$other_model->analysis_project, 'unassigned model not assigned to');

my $wrong_class_profile_item = add_config($other_analysis_project, 'Genome::Model::SomaticValidation');
my $cmd6 = $class->create(
    profile_item => $wrong_class_profile_item,
    models => [$other_model],
);
isa_ok($cmd6, $class, 'created command');
ok(!$cmd6->execute, 'command fails when attempting to assign to a config profile item for the wrong model type');
ok(!$other_model->analysis_project, 'unassigned model not assigned to');

$other_profile_item->status('disabled');
my $cmd7 = $class->create(
    profile_item => $other_profile_item,
    models => [$other_model],
);
isa_ok($cmd7, $class, 'created command');
ok($cmd7->execute, 'command succeeds in "normal" case');
is($other_model->analysis_project, $other_analysis_project, 'model assigned appropriately');
is($other_model->config_profile_items, $other_profile_item, 'model linked to correct configuration');


sub add_config {
    my $analysis_project = shift;
    my $model_class = shift;

    my $config = <<EOFILE
models:
    "$model_class":
        processing_profile_id: 1
EOFILE
;
    my $config_file = Genome::Sys->create_temp_file_path;
    Genome::Sys->write_file($config_file, $config);

    my $add_cmd = Genome::Config::AnalysisProject::Command::AddConfigFile->create(
        analysis_project => $analysis_project,
        config_file => $config_file,
        store_only => 1,
        reprocess_existing => 0,
    );
    my $profile_item = $add_cmd->execute();
    isa_ok($profile_item, 'Genome::Config::Profile::Item', 'created config');
    return $profile_item;
}
