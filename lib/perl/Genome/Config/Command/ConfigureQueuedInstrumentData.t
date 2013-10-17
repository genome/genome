#!/usr/bin/env genome-perl
use strict;
use warnings;

$ENV{UR_DBI_NO_COMMIT} = 1;

use Test::More;
use above "Genome";
use Carp::Always;

use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::AnalysisProject;
use Genome::Test::Factory::Individual;
use Genome::Test::Factory::Library;
use Genome::Test::Factory::Sample;
use Genome::Test::Factory::Taxon;

no warnings qw(redefine);
use UR::Context::Transaction;
*UR::Context::Transaction::commit = sub { return 1; };
use warnings;

my $class = 'Genome::Config::Command::ConfigureQueuedInstrumentData';
use_ok($class);

#new model
my ($rna_instrument_data, $model_types) = generate_rna_seq_instrument_data();
build_and_run_cmd($rna_instrument_data);
assert_succeeded($rna_instrument_data, $model_types);
is($rna_instrument_data->models->auto_assign_inst_data, 1, 'it should default to setting auto assign inst data to 1');

#existing model
($rna_instrument_data, $model_types) = generate_rna_seq_instrument_data();
my $config_hash = _rna_seq_config_hash();
delete $config_hash->{instrument_data_properties};
$config_hash->{subject} = $rna_instrument_data->sample;
$config_hash->{target_region_set_name} = $rna_instrument_data->target_region_set_name;
$config_hash->{auto_assign_inst_data} = 1;
my $rna_model = Genome::Model::RnaSeq->create(%{$config_hash});
build_and_run_cmd($rna_instrument_data);
assert_succeeded($rna_instrument_data, $model_types);
is($rna_model->instrument_data->id, $rna_instrument_data->id, 'it assigned the instrument data to the existing model');

#existing model without auto assign inst data set
($rna_instrument_data, $model_types) = generate_rna_seq_instrument_data();
my $config_hash_no_auto_assign = _rna_seq_config_hash();
delete $config_hash_no_auto_assign->{instrument_data_properties};
$config_hash_no_auto_assign->{subject} = $rna_instrument_data->sample;
$config_hash_no_auto_assign->{target_region_set_name} = $rna_instrument_data->target_region_set_name;
$config_hash_no_auto_assign->{auto_assign_inst_data} = 0;
my $rna_model_no_auto_assign = Genome::Model::RnaSeq->create(%{$config_hash_no_auto_assign});
build_and_run_cmd($rna_instrument_data);
assert_succeeded($rna_instrument_data, $model_types);
is($rna_model_no_auto_assign->instrument_data, undef , 'it did not assign the instrument data to the existing model');

#model with pairing
my ($data1, $data2, $model_types_somval, $analysis_project) = _generate_som_val_instrument_data();
Genome::Config::AnalysisProject::SubjectPairing->create(
    analysis_project => $analysis_project,
    control_subject => @$data1[1],
    experimental_subject => @$data2[1]
);
Genome::Config::AnalysisProject::SubjectPairing->create(
    analysis_project => $analysis_project,
    control_subject => @$data2[1],
    experimental_subject => @$data1[1]
);

build_and_run_cmd(@$data1[0], @$data2[0]);
assert_succeeded(@$data1[0], $model_types_somval);
assert_succeeded(@$data2[0], $model_types_somval);
my @models = $analysis_project->models;
ok(@models, 'it registers created models with the analysis_project');
ok(scalar(@models) == 2, 'it creates one model per SubjectPairing');

#model with pairing - but none set
($data1, $data2, $model_types_somval, $analysis_project) = _generate_som_val_instrument_data();
build_and_run_cmd(@$data1[0], @$data2[0]);
assert_failed(@$data1[0], 'Found no pairing information');
assert_failed(@$data2[0], 'Found no pairing information');

#inst data with no ap
my $inst_data_without_a_project = Genome::Test::Factory::InstrumentData::Solexa->setup_object();
build_and_run_cmd($inst_data_without_a_project);
assert_failed($inst_data_without_a_project, "doesn't have an analysis project!");
done_testing();


sub assert_failed {
    my $inst_data = shift;
    my $error = shift;
    ok($inst_data->tgi_lims_status eq 'failed', 'it should mark the inst data as failed');
    ok($inst_data->attributes(attribute_label => 'tgi_lims_fail_count'), 'it should set the fail count');
    ok($inst_data->attributes(attribute_label => 'tgi_lims_fail_message'), 'it should set the fail message');
    if ($error) {
        ok($inst_data->attributes(attribute_label => 'tgi_lims_fail_message')->value =~ $error, 'it should set the correct error message');
    }
}

sub assert_succeeded {
    my $inst_data = shift;
    my $model_types = shift;

    ok($inst_data->tgi_lims_status eq 'processed', 'it should mark the inst data as succeeded');
    ok(!$inst_data->attributes(attribute_label => 'tgi_lims_fail_count'), 'it should remove the fail count');
    for my $model_instance ($inst_data->models) {
        ok($model_instance->build_requested, 'it sets build requested on constructed models');
    }
    for my $model_type (@$model_types) {
        ok(
            scalar(grep { $_->class eq $model_type } $inst_data->models),
            'it creates a model of the expected type and associates it with the inst data'
        );
    }
}

sub build_and_run_cmd {
    my @inst_data = @_;
    my $cmd = Genome::Config::Command::ConfigureQueuedInstrumentData->create(
        instrument_data => [@inst_data],
    );

    ok($cmd, 'created the command successfully');
    ok($cmd->execute(), 'command ran successfully');
    return;
}

sub generate_rna_seq_instrument_data {
    my $inst_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
        fwd_read_length => 20,
        fwd_clusters => 10,
        rev_read_length => 20,
        rev_clusters => 10,
        target_region_set_name => 'turkey',
        run_type => 'Paired',
    );

    my $ap = Genome::Test::Factory::AnalysisProject->setup_object(
        config_hash => {
            'Genome::Model::RnaSeq' => _rna_seq_config_hash()
        }
    );

    $inst_data->analysis_project_id($ap->id);
    return ($inst_data, ['Genome::Model::RnaSeq']);
}


sub _rna_seq_config_hash {
    return {
        processing_profile_id       => 2819506,
        annotation_build_id         => 124434505,
        reference_sequence_build_id => 106942997,
        instrument_data_properties  => {
            subject => 'sample',
            target_region_set_name => 'target_region_set_name'
        }
    };
}

sub _generate_som_val_instrument_data {
    my $sample_1 = Genome::Test::Factory::Sample->setup_object(
        extraction_type => 'genomic_dna',
    );

    my $sample_2 = Genome::Test::Factory::Sample->setup_object(
        extraction_type => 'genomic_dna',
        source_id => $sample_1->source->id,
    );

    my $library_1 = Genome::Test::Factory::Library->setup_object(
        sample_id => $sample_1->id,
    );

    my $library_2 = Genome::Test::Factory::Library->setup_object(
        sample_id => $sample_2->id,
    );

    my $inst_data_1 = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
        library_id => $library_1->id,
        fwd_read_length => 20,
        fwd_clusters => 10,
        rev_read_length => 20,
        rev_clusters => 10,
        run_type => 'Paired',
    );

    my $inst_data_2 = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
        library_id => $library_2->id,
        fwd_read_length => 20,
        fwd_clusters => 10,
        rev_read_length => 20,
        rev_clusters => 10,
        run_type => 'Paired',
    );

    my $ap = Genome::Test::Factory::AnalysisProject->setup_object(
        config_hash => {
            'Genome::Model::SomaticValidation' => _som_val_config_hash()
        }
    );

    $inst_data_1->analysis_project_id($ap->id);
    $inst_data_2->analysis_project_id($ap->id);
    return ([$inst_data_1, $sample_1], [$inst_data_2, $sample_2], ['Genome::Model::SomaticValidation'], $ap);
}

sub _som_val_config_hash {
    return {
        processing_profile_id       => 2656116,
        annotation_build_id         => 124434505,
        reference_sequence_build_id => 106942997,
    };
}
