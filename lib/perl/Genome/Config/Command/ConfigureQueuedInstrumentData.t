#!/usr/bin/env genome-perl
use strict;
use warnings;

$ENV{UR_DBI_NO_COMMIT} = 1;

use Test::More;
use above "Genome";
use Carp::Always;

use Genome::TestObjGenerator::AnalysisProject;
use Genome::TestObjGenerator::InstrumentData::Solexa;

my $class = 'Genome::Config::Command::ConfigureQueuedInstrumentData';
use_ok($class);

##construct refalign
#my $ref_align_instrument_data = ref_align_instrument_data();
#build_and_run_cmd($ref_align_instrument_data);
#assert_succeeded($ref_align_instrument_data);

##existing refalign


##existing model
##construct rnaseq
##construct somval - multiple subject pairings

##inst data with no config

##inst data with no ap

#done_testing();

#sub assert_failed {
    #my $inst_data = shift;
    #ok($inst_data->tgi_lims_status eq 'failed', 'it should mark the inst data as failed');
    #ok($inst_data->attributes(attribute_label => 'tgi_lims_fail_count'), 'it should set the fail count');
    #ok($inst_data->attributes(attribute_label => 'tgi_lims_fail_message'), 'it should set the fail message');
#}

#sub assert_succeeded {
    #my $inst_data = shift;
    #ok($inst_data->tgi_lims_status eq 'processed', 'it should mark the inst data as failed');
    #ok(!$inst_data->attributes(attribute_label => 'tgi_lims_fail_count'), 'it should remove the fail count');
#}

#sub build_and_run_cmd {
    #my $inst_data = shift;
    #my $cmd = Genome::Config::Command::ConfigureQueuedInstrumentData->create(
        #instrument_data => $inst_data,
    #);

    #ok($cmd, 'created the command successfully');
    #ok($cmd->execute(), 'command ran successfully');
    #return;
#}

#sub ref_align_instrument_data {
#}

#sub rna_seq_instrument_data {
    #my $inst_data = Genome::TestObjGenerator::InstrumentData::Solexa->setup_object();
    #my $ap = Genome::TestObjGenerator::AnalysisProject->setup_object(
        #config_hash => {
            #'Genome::Model::RnaSeq' => {
                #processing_profile_id       => 2819506,
                #annotation_build_id         => 124434505,
                #reference_sequence_build_id => 106942997,
                #instrument_data_properties  => {
                    #subject => 'sample',
                    #target_region_set_name => 'target_region_set_name'
                #},
            #}
        #}
    #);

    #$inst_data->analysis_project_id($ap->id);
    #return $inst_data;
#}

#sub somatic_validation_instrument_data {
#}

#sub instrument_data_with_no_config {
#}

#sub instrument_data_with_no_analysis_project {
#}


