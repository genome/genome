#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $SIG{__DIE__} = sub{ Carp::confess(@_); };
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
    $ENV{UR_COMMAND_DUMP_DEBUG_MESSAGES} = 1;
};

use above 'Genome';
use Genome::Test::Factory::DiskAllocation;
use Genome::Test::Factory::InstrumentData::Solexa;

use Test::Exception;
use Test::More tests => 4;

my $class = 'Genome::Config::AnalysisProject::Command::ReplaceModel';
use_ok($class) or die;

my $model_class = 'Genome::Model::Tester';
my $analysis_project = Genome::Config::AnalysisProject->create(
    run_as => Genome::Sys->username,
    status => 'Completed',
);
class Genome::Model::Tester { is => 'Genome::Model', };
my ($model);

subtest "setup" => sub{
    plan tests => 9;

    # create a pp and profile item
    for my $id ( -11, -22 ) {
        my $pp = Genome::ProcessingProfile::Tester->create(
            id => $id,
            name => $id,
        );
        ok($pp, 'create pp');
        my $profile_item = Genome::Config::Profile::Item->create(
            id => $id,
            analysis_project => $analysis_project,
            status => 'active',
        );
        ok($profile_item, 'create profile item');
        my $profile_alloc = Genome::Test::Factory::DiskAllocation->generate_obj(
            owner => $profile_item,
            mount_path => Genome::Sys->create_temp_file_path,
        );
    my $data = <<EOFILE
rules:
  sequencing_platform: solexa

models:
  "$model_class":
    processing_profile_id: $id
    instrument_data_properties:
      subject: sample
EOFILE
;
        Genome::Sys->write_file( File::Spec->join($profile_alloc->absolute_path, 'config.yml'), $data );
        ok(-s $profile_item->file_path, 'created profile item file path');
    }

    my $instdata = Genome::Test::Factory::InstrumentData::Solexa->setup_object;
    $model = $model_class->create(
        run_as => Genome::Sys->username,
        subject => $instdata->sample,
        processing_profile_id => -11,
        auto_assign_inst_data => 1
    );
    $model->add_instrument_data($instdata);
    ok($model->instrument_data, 'added model instrument data');

    my $profile_item =  Genome::Config::Profile::Item->get(-11);
    $analysis_project->add_model_bridge(
        model => $model, config_profile_item => $profile_item,
    );
    is($model->analysis_project, $analysis_project, 'model assigned to anp');
    is($model->config_profile_item, $profile_item, 'model has config item');

};

subtest "fails" => sub{
    plan tests => 3;

    my $profile_item = Genome::Config::Profile::Item->get(-11);

    throws_ok(
        sub{
            $class->execute(
                new_profile_item => $profile_item,
                model => $model,
            ); },
        qr/status \(active\) must be 'inactive'/,
        'fails when profile item is active',
    );

    $profile_item->status('inactive');
    throws_ok(
        sub{
            $class->execute(
                new_profile_item => $profile_item,
                model => $model,
            ); },
        qr/status \(Completed\) is not a 'current'/,
        'fails when anp is not current',
    );

    $analysis_project->status("In Progress");
    throws_ok(
        sub{
            $class->execute(
                new_profile_item => $profile_item,
                model => $model,
            ); },
        qr/No model found or created/,
        'fails when no model found/created',
    );

};

subtest "replace model" => sub{
    plan tests => 5;

    $analysis_project->status("In Progress");
    my $new_profile_item = Genome::Config::Profile::Item->get(-22);
    $new_profile_item->status('inactive');

    my $cmd = $class->create(
        new_profile_item => $new_profile_item,
        model => $model,
    );

    lives_ok( sub{ $cmd->execute }, 'execute lives');
    ok($cmd->result, 'result is successful');
    ok($cmd->new_model, 'created new model');
    is($cmd->new_model->analysis_project, $analysis_project, 'new model assigned to anp');
    is($cmd->new_model->config_profile_item, $new_profile_item, 'new model has new config item');

};

done_testing();
