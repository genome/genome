#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';
use Test::More tests => 4;

use Genome::Test::Factory::AnalysisProject;
use Genome::Test::Factory::Model::SomaticValidation;
use Genome::Test::Factory::Build;
use Genome::Test::Factory::Sample;
use Genome::Utility::Test;
Genome::Report::Email->silent();

my $pkg = 'Genome::Model::SomaticValidation::Command::ValidateSvs::CreateAssembledContigReference';
use_ok($pkg) or die();

my $test_dir = Genome::Utility::Test->data_dir_ok($pkg);

my $anp = Genome::Test::Factory::AnalysisProject->setup_object();

my $env_file = Genome::Sys->create_temp_file_path;
Genome::Sys->write_file($env_file, 'not_a_real_config_file: 1');
Genome::Config::AnalysisProject::Command::AddEnvironmentFile->execute(environment_file => $env_file, analysis_project => $anp);

my $model = Genome::Test::Factory::Model::SomaticValidation->setup_object();
$anp->add_model_bridge(model_id => $model->id);
my $build = Genome::Test::Factory::Build->setup_object(model_id => $model->id, status => 'Running');
my $sample = Genome::Test::Factory::Sample->setup_object(source_id => $model->subject_id);
$build->tumor_sample($sample);
my $cmd = $pkg->create(
    build_id => $build->id,
    skip => 0,
);
isa_ok($cmd, $pkg, 'created command');

my $output_dir = $cmd->output_dir;
Genome::Sys->create_directory($output_dir);
Genome::Sys->symlink_directory($test_dir, $output_dir);

ok($cmd->execute, 'command executes');


