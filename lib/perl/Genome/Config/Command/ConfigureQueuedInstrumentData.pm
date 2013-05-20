package Genome::Config::Command::ConfigureQueuedInstrumentData;

use strict;
use warnings;

class Genome::Config::Command::ConfigureQueuedInstrumentData {
    is => 'Command::V2'
};

sub help_brief {
    return 'Assign instrument data with an analysis project to models';
}

sub help_synopsis {
    return <<'EOS'
EOS
}

sub help_detail {
    return <<EOS
EOS
}

sub execute {
    my $self = shift;

    $self->_lock();

    my @instrument_data = $self->_get_instrument_data_to_process() || return 1;

    for my $current_inst_data (@instrument_data) {
        my $config = $self->_get_configuration_object_for_instrument_data($current_inst_data);
        my $hashes = $self->_prepare_configuration_hashes_for_instrument_data($current_inst_data, $config);
        for(values %$hashes) {
            my $model = $self->_get_model_for_config_hash($_);
            $model->add_instrument_data($current_inst_data);
            $model->build_requested(1);
        }
    }

    UR::Context->commit();
}

sub _get_model_for_config_hash {
    #hashref
    my $config = shift;

    my $m = Genome::Model->get(%$config);
    return $m || Genome::Model->create(%$config);
}

#SAMPLE FORMAT of YAML files:

#"Genome::Model::ReferenceAlignment":
#    processing_profile_id: 123451
#    annotation_reference_build_id: 12314
#    dbsnp_build_id: 123124
#    reference_sequence_build_id: 12314
#"Genome::Model::RnaSeq":
#    processing_profile_id: 123451
#    annotation_build_id: 12314
#    reference_sequence_build_id: 12314
#this gets read in as a hash(ref) of hash(refs) - one hash per model type
sub _prepare_configuration_hashes_for_instrument_data {
    #given a piece of ID and a configuration object, prepare a hash
    #to pass to model get or create
    my ($inst_data, $config_obj) = @_;
    #eventually this will need to support multiple references
    my $config_hash = $config_obj->get_config(
        taxon => $inst_data->species_name,
        type => $inst_data->extraction_type,
    );
    for(keys %$config_hash) {
        $config_hash->{$_}->{subclass_name} = $_;
        $config_hash->{$_}->{subject} = $inst_data->subject;
        $config_hash->{$_}->{target_region_set_name} = $inst_data->target_region_set_name;
        #not sure how to derive this yet
        #$config_hash->{$_}->{genotype_microarray} = 
    }
    return $config_hash;
}

sub _get_instrument_data_to_process {
    my $self = shift;

    #TODO switch to a boolexpr maybe? there's got to be a better way to do this
    return grep { $_->analysis_project_id } Genome::InstrumentData->get(
        'tgi_lims_status' => [qw/ new failed /],
        -hint => [ 'sample', 'sample.source', 'sample.source.taxon', ],
    );
}

sub _get_configuration_object_for_instrument_data {
    my $self = shift;
    my $instrument_data = shift;
    die('You must provide a single piece of instrument data!') unless $instrument_data;

    my $analysis_project_id = $instrument_data->analysis_project_id;
    my $analysis_project = Genome::Config::AnalysisProject->get($analysis_project_id)
        || die("$analysis_project_id doesn't seem to be a valid Genone::Config::AnalysisProject ID!");
    return $analysis_project->get_configuration_reader();
}

sub _lock {
    my $lock_var = $ENV{GENOME_LOCK_DIR} . '/genome_config_command_configure-queued-instrument-data/lock';
    my $lock = Genome::Sys->lock_resource(resource_lock => $lock_var, max_try => 1);

    die('Unable to acquire the lock! Is ConfigureQueuedInstrumentData already running or did it exit uncleanly?')
        unless $lock;

    UR::Context->current->add_observer(
        aspect => 'commit',
        callback => sub{
            Genome::Sys->unlock_resource(resource_lock=>$lock);
        }
    );
}

1;
