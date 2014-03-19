#!/usr/bin/env genome-perl
use strict;
use warnings;

$ENV{UR_DBI_NO_COMMIT} = 1;

use Test::More;
use above "Genome";
use Carp::Always;

use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::InstrumentData::Imported;
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
is_deeply([$rna_model->instrument_data], [], 'it does not assign the instrument data to an existing model outside the analysis project');

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
my $subject_mapping = Genome::Config::AnalysisProject::SubjectMapping->create(
    analysis_project => $analysis_project
);

Genome::Config::AnalysisProject::SubjectMapping::Subject->create(
    label => 'control_subject',
    subject_mapping => $subject_mapping,
    subject => @$data1[1],
);
Genome::Config::AnalysisProject::SubjectMapping::Subject->create(
    label => 'subject',
    subject_mapping => $subject_mapping,
    subject => @$data2[1]->source,
);
Genome::Config::AnalysisProject::SubjectMapping::Subject->create(
    label => 'experimental_subject',
    subject_mapping => $subject_mapping,
    subject => @$data2[1],
);

my $subject_mapping2 = Genome::Config::AnalysisProject::SubjectMapping->create(
    analysis_project => $analysis_project
);
Genome::Config::AnalysisProject::SubjectMapping::Subject->create(
    label => 'control_subject',
    subject_mapping => $subject_mapping2,
    subject => @$data2[1],
);
Genome::Config::AnalysisProject::SubjectMapping::Subject->create(
    label => 'subject',
    subject_mapping => $subject_mapping2,
    subject => @$data1[1]->source,
);
Genome::Config::AnalysisProject::SubjectMapping::Subject->create(
    label => 'experimental_subject',
    subject_mapping => $subject_mapping2,
    subject => @$data1[1],
);

build_and_run_cmd(@$data1[0], @$data2[0]);
assert_succeeded(@$data1[0], $model_types_somval);
assert_succeeded(@$data2[0], $model_types_somval);
my @models = $analysis_project->models;
ok(@models, 'it registers created models with the analysis_project');
ok(scalar(@models) == 2, 'it creates one model per SubjectMapping');

#model with pairing - but none set
($data1, $data2, $model_types_somval, $analysis_project) = _generate_som_val_instrument_data();
build_and_run_cmd(@$data1[0], @$data2[0]);
assert_failed(@$data1[0], 'Found no mapping information');
assert_failed(@$data2[0], 'Found no mapping information');

#inst data with no ap
my $inst_data_without_a_project = Genome::Test::Factory::InstrumentData::Solexa->setup_object();
my $cmd = build_and_run_cmd($inst_data_without_a_project);
ok($cmd->status_message =~ /Found no items to process/, 'no analysis project is a no-op');

#skipped data
($rna_instrument_data, $model_types) = generate_rna_seq_instrument_data();
$rna_instrument_data->ignored(1);
build_and_run_cmd($rna_instrument_data);
assert_skipped($rna_instrument_data);

#lane_qc and genotype handling
my ($sans_data, $plus_data, $model_types) = _generate_lane_qc_instrument_data();
build_and_run_cmd($sans_data, $plus_data);
assert_succeeded($plus_data, $model_types);
my ($plus_model) = $plus_data->models;
ok($plus_model->build_requested, 'Lane QC model with genotype data has build requested');
ok($plus_model->genotype_microarray, 'Lane QC model with genotype data has genotype_microarray model');
my ($sans_model) = $sans_data->models;
ok(!$sans_model->build_requested, 'Lane QC model without genotype data does not have build requested');
ok(!$sans_model->genotype_microarray, 'Lane QC model without genotype data does not have genotype_microarray model');

done_testing();

sub assert_skipped {
    my $inst_data = shift;

    my ($bridge) = $inst_data->analysis_project_bridges;
    ok($bridge->status eq 'skipped', 'it should mark the inst data as skipped');
    ok($bridge->reason, 'a reason for being skipped was specified');
    is($bridge->fail_count, 0, 'a fail count was not set');
}

sub assert_failed {
    my $inst_data = shift;
    my $error = shift;
    my ($bridge) = $inst_data->analysis_project_bridges;
    ok($bridge->status eq 'failed', 'it should mark the inst data as failed');
    ok($bridge->fail_count, 'it should set the fail count');
    ok($bridge->reason, 'it should set the fail message');
    if ($error) {
        ok($bridge->reason =~ $error, 'it should set the correct error message');
    }
}

sub assert_succeeded {
    my $inst_data = shift;
    my $model_types = shift;
    my ($bridge) = $inst_data->analysis_project_bridges;

    ok($bridge->status eq 'processed', 'it should mark the inst data as succeeded');
    is($bridge->fail_count, 0, 'it should remove the fail count');
    for my $model_instance ($inst_data->models) {
        ok($model_instance->build_requested, 'it sets build requested on constructed models');
        is($model_instance->user_name, 'apipe-builder');
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
    return $cmd;
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

    Genome::Config::AnalysisProject::InstrumentDataBridge->create(
        instrument_data => $inst_data,
        analysis_project => $ap,
    );

    return ($inst_data, ['Genome::Model::RnaSeq']);
}


sub _rna_seq_config_hash {
    return {
        processing_profile_id       => 2819506,
        annotation_build_id         => 124434505,
        reference_sequence_build_id => 106942997,
        user_name => 'apipe-builder',
        instrument_data_properties  => {
            subject => 'sample',
            target_region_set_name => 'target_region_set_name'
        }
    };
}

my $som_val_project;
sub _generate_som_val_instrument_data {
    $som_val_project ||= Genome::Project->create(name => sprintf('test %s %s',  __FILE__));

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

    map $som_val_project->add_part(role => 'instrument_data', entity => $_), ($inst_data_1, $inst_data_2);

    my $ap = Genome::Test::Factory::AnalysisProject->setup_object(
        config_hash => {
            'Genome::Model::SomaticValidation' => _som_val_config_hash()
        }
    );

    for my $inst_data ($inst_data_1, $inst_data_2) {
        Genome::Config::AnalysisProject::InstrumentDataBridge->create(
            instrument_data => $inst_data,
            analysis_project => $ap,
        );
    }

    return ([$inst_data_1, $sample_1], [$inst_data_2, $sample_2], ['Genome::Model::SomaticValidation'], $ap);
}

sub _som_val_config_hash {
    return {
        processing_profile_id       => 2656116,
        annotation_build_id         => 124434505,
        reference_sequence_build_id => 106942997,
        user_name                   => 'apipe-builder',
    };
}

sub _generate_lane_qc_instrument_data {
    my $genotype_sample = Genome::Test::Factory::Sample->setup_object();

    my $genotype_library => Genome::Test::Factory::Library->setup_object(
        sample_id => $genotype_sample, 
    );

    my $genotype_data = Genome::Test::Factory::InstrumentData::Imported->setup_object(
        library => $genotype_library,
    );

    #Lane QC stuff sans genotype data
    my $sans_sample = Genome::Test::Factory::Sample->setup_object();

    my $sans_library = Genome::Test::Factory::Library->setup_object(
        sample_id => $sans_sample, 
    );

    my $sans_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
        library => $sans_library,
        run_name => 'sans',
        subset_name => 3,
    );

    #Lane QC stuff plus genotype data
    my $plus_sample = Genome::Test::Factory::Sample->setup_object(
        extraction_type => 'genomic_dna',
        default_genotype_data_id => , $genotype_data->id,
    );

    my $plus_library = Genome::Test::Factory::Library->setup_object(
        sample_id => $plus_sample, 
    );

    my $plus_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
        library => $plus_library,
        run_name => 'plus',
        subset_name => 2,
    );

    my $genotype_microarray_model = Genome::Model::GenotypeMicroarray->create(
        processing_profile_id => '2166945',
        dbsnp_build_id => 127786607,
        subject => $genotype_sample,
        instrument_data => [$genotype_data],
    );

    my $tmp_dir = Genome::Sys->create_temp_directory;
    my $gmb = Genome::Model::Build->create(model_id => $genotype_microarray_model->id, data_directory => $tmp_dir);
    $gmb->success();


    for my $inst_data ($sans_data, $plus_data) {
        my $ap = Genome::Test::Factory::AnalysisProject->setup_object(
            config_hash => {
                'Genome::Model::ReferenceAlignment' => _ref_align_config_hash()
            }
        );

        Genome::Config::AnalysisProject::InstrumentDataBridge->create(
            instrument_data => $inst_data,
            analysis_project => $ap,
        );
    }

    return ($sans_data, $plus_data, ['Genome::Model::ReferenceAlignment']); 
}


sub _ref_align_config_hash {
    return {
        processing_profile_id => '2653572',
        reference_sequence_build_id => '106942997',
        user_name                   => 'apipe-builder',
        auto_assign_inst_data => 0,
        instrument_data_properties  => {
            subject => 'sample',
            instrument_data => ['__self__'],
        }
    };
}
