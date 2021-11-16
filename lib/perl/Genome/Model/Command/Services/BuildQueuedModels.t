#!/usr/bin/env genome-perl

use strict;
use warnings;

use Sub::Install qw();
use Test::More tests => 14;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Genome::Test::Factory::AnalysisProject;
use Genome::Test::Factory::Model::CwlPipeline;

my $class = 'Genome::Model::Command::Services::BuildQueuedModels';
use_ok($class);

setup_overrides();

my @models = setup_data();
is(_build_requested_count(@models), 12, 'setup test data');


my $cmd = $class->create(
    scheduled_max => 3,
    per_analysis_project_scheduled_max => 1,
    running_max => 1000000,
    start_limit => 1000000,
    submit_jobs => 0,
);
isa_ok($cmd, $class, 'created command');
ok($cmd->execute, 'executed command');
is(_build_requested_count(@models), 9, 'only scheduled up to scheduled max');
my @build = Genome::Model::Build::CwlPipeline->get(status => 'Scheduled', run_by => Genome::Sys->username);
is(scalar(@build), 3, 'scheduled max number of builds');
my %anp;
for (@build ) {
    $anp{ $_->model->analysis_project->id }++;
}
is(scalar keys(%anp), 3, 'started builds across three AnPs');

#so UR will have correct behavior on future iterators
UR::Context->commit();
Genome::Model->unload();
Genome::Model::Set->unload();

my $cmd_no_new = $class->create(
    scheduled_max => 2,
    per_analysis_project_scheduled_max => 1000000,
    running_max => 1000000,
    start_limit => 1000000,
    submit_jobs => 0,
);
isa_ok($cmd_no_new, $class);
ok($cmd_no_new->execute, 'executed command');
is(_build_requested_count(@models), 9, 'build requested count unchanged');

#so UR will have correct behavior on future iterators
UR::Context->commit();
Genome::Model->unload();
Genome::Model::Set->unload();

my $cmd_start_limit = $class->create(
    scheduled_max => 1000000,
    per_analysis_project_scheduled_max => 1000000,
    running_max => 1000000,
    start_limit => 5,
    submit_jobs => 0,
);
isa_ok($cmd_start_limit, $class, 'created command');
ok($cmd_start_limit->execute, 'executed command');
is(_build_requested_count(@models), 4, 'only scheduled up to scheduled max');
 @build = Genome::Model::Build::CwlPipeline->get(status => 'Scheduled', run_by => Genome::Sys->username);
is(scalar(@build), 8, 'scheduled max number of builds');

done_testing();

sub setup_overrides {
    Genome::Sys->class;
    my $name = 'fake_user_for_testing-BQM' . time();
    Sub::Install::reinstall_sub({
        into => 'Genome::Sys',
        as => 'username',
        code => sub { return $name },
    });

    Genome::Model::Build->class;
    Sub::Install::reinstall_sub({
        into => 'Genome::Model::Build',
        as => '_launch',
        code => sub { return 1; },
    });
}

sub setup_data {
    my @models;
    for (1..4) {
        my $anp = Genome::Test::Factory::AnalysisProject->setup_object;
        my $env_file = Genome::Sys->create_temp_file_path;
        Genome::Sys->write_file($env_file, 'not_a_real_config_file: 1');
        Genome::Config::AnalysisProject::Command::AddEnvironmentFile->execute(environment_file => $env_file, analysis_project => $anp);

        for (1..5) {
            my $m = Genome::Test::Factory::Model::CwlPipeline->setup_object(
                run_as => Genome::Sys->username,
            );
            $anp->add_model_bridge(model_id => $m->id);
            $m->build_requested($_ % 2);

            push @models, $m;
        }
    }

    #so UR will have correct behavior on future iterators
    UR::Context->commit();
    Genome::Model->unload();

    return @models;
}

sub _build_requested_count {
    return scalar grep { $_->build_requested } @_;
}

