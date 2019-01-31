#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::Exception;
use Test::More tests => 57;
use Test::Deep qw(cmp_bag);


use above 'Genome';

use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::Library;
use Genome::Test::Factory::ProcessingProfile::CwlPipeline;
use Genome::Test::Factory::Sample;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

my $analysis_project = _setup_analysis_project();
my $instrument_data = _setup_instrument_data($analysis_project);

my $profile = Genome::Config::Profile->create_from_analysis_project($analysis_project);
isa_ok($profile, 'Genome::Config::Profile', 'created profile');

my %models;

dies_ok { $profile->process_models_for_instrument_data($instrument_data->[0][0]); } 'fails when no mapping defined';

my $subject_mapping = Genome::Config::AnalysisProject::SubjectMapping->create(
    analysis_project => $analysis_project
);

my $sample = $instrument_data->[0][0]->sample;

$subject_mapping->add_subject_bridge(subject => $sample, label => 'best_sample');
$subject_mapping->add_input(key => 'some_key', value => 'some_value');


for my $set (@$instrument_data) {
    my @next_models = $profile->process_models_for_instrument_data($set->[0]);
    is(scalar(@next_models), 2, 'resulted in one model as expected');
    ok(!exists $models{$next_models[0]->id}, 'model was newly created');
    $models{$next_models[0]->id} = $next_models[0];
    $models{$next_models[1]->id} = $next_models[1];

    my @next_models2 = $profile->process_models_for_instrument_data($set->[1]);
    is(scalar(@next_models), 2, 'resulted in one model as expected');
    ok(exists $models{$next_models2[0]->id}, 'used existing model');
    ok(exists $models{$next_models2[1]->id}, 'used existing model');
    $models{$next_models2[0]->id} = $next_models2[0];
    $models{$next_models2[1]->id} = $next_models2[1];

    my @i = $next_models2[0]->instrument_data;
    cmp_bag(\@i, $set, 'model has expected instrument data assigned');
    my @i2 = $next_models2[1]->instrument_data;
    cmp_bag(\@i2, $set, 'model has expected instrument data assigned');
}

is(scalar keys %models, 6, 'created expected number of models');

for my $m (values %models) {
    my @inputs = $m->inputs;
    my %inputs = map { $_->name => $_->value_id } @inputs;

    ok(exists $inputs{this_uses_mappings}, 'test mapping indicator set');

    if ($inputs{this_uses_mappings}) {
        is($inputs{best_sample}, $sample->id, 'sample attached');
        is($inputs{some_key}, 'some_value', 'extra key attached');
    } else {
        #test that static inputs were included
        is($inputs{bird1}, 'turkey', 'set first bird input');
        is($inputs{bird2}, 'eagle', 'set second bird input');
        is($inputs{bird3}, 'dove', 'set third bird input');
    }

    #test that dynamic inputs were included
    ok($inputs{library}, 'library input is set');
    ok($inputs{flow_cell}, 'flow_cell input is set');
}


sub _setup_instrument_data {
    my $anp = shift;

    my $sample = Genome::Test::Factory::Sample->setup_object();
    my $lib1 = Genome::Test::Factory::Library->setup_object(sample => $sample);
    my $lib2 = Genome::Test::Factory::Library->setup_object(sample => $sample);

    my @i_lib1 = map { Genome::Test::Factory::InstrumentData::Solexa->setup_object(library => $lib1, flow_cell_id => 'TESTCCXX123456789', lane => $_) } (1..2);
    my @i_lib2 = map { Genome::Test::Factory::InstrumentData::Solexa->setup_object(library => $lib2, flow_cell_id => 'TESTCCXX123456789', lane => $_) } (3..4);
    my @i_lib2_otherflowcell = map { Genome::Test::Factory::InstrumentData::Solexa->setup_object(library => $lib1, flow_cell_id => 'TESTAAXX987654321', lane => $_) } (5..6);

    for (@i_lib1, @i_lib2, @i_lib2_otherflowcell) {
        Genome::Config::AnalysisProject::InstrumentDataBridge->create(
            instrument_data_id => $_->id,
            analysis_project_id => $anp->id,
        );
    }

    return [\@i_lib1, \@i_lib2, \@i_lib2_otherflowcell];
}

sub _setup_analysis_project {
    my $anp = Genome::Config::AnalysisProject->create(
        name => 'Test Project for input data test',
        run_as => 'nobody'
    );

    my $pp = Genome::Test::Factory::ProcessingProfile::CwlPipeline->setup_object;

    my $config_file = _config_file($pp);
    my $menu_item = Genome::Config::AnalysisMenu::Item->create(
        id => '-12345678',
        file_path => $config_file,
        name => 'input data test menu item',
    );

    my $item = Genome::Config::Profile::Item->create(
        analysis_project => $anp,
        analysis_menu_item => $menu_item,
        status => 'active',
    );

    return $anp;
}


sub _config_file {
    my $pp = shift;
    my $ppid = $pp->id;

    my $contents = <<YML
rules:
  sequencing_platform: solexa
models:
  'Genome::Model::CwlPipeline':
     - processing_profile_id: $ppid
       input_data:
           this_uses_mappings: 0
           bird1: turkey
           bird2: eagle
           bird3: dove
       instrument_data_properties:
           input_data:
               library: library
               flow_cell: flow_cell_id
           subject: sample
     - processing_profile_id: $ppid
       input_data_requires_subject_mapping: 1
       input_data:
           this_uses_mappings: 1
       instrument_data_properties:
           subject: sample
           input_data:
               library: library
               flow_cell: flow_cell_id
YML
    ;

    my $file_path = join('', Genome::Sys->create_temp_file_path(), '.yml');

    Genome::Sys->write_file($file_path, $contents);

    return $file_path;
}
