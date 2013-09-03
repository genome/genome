package Genome::Config::Command::ConfigureQueuedInstrumentData;

use strict;
use warnings;

class Genome::Config::Command::ConfigureQueuedInstrumentData {
    is => 'Command::V2',
    has_optional => [
        instrument_data => {
            is          => 'Genome::InstrumentData',
            is_many     => 1,
            doc         => '[Re]process these instrument data.',
        },
    ],
};

sub help_brief {
    return 'Assign instrument data with an analysis project to models';
}

sub help_synopsis {
    return <<EOS
EOS
}

sub help_detail {
    return <<EOS
EOS
}

sub execute {
    my $self = shift;

    $self->_lock();

    my @instrument_data = $self->_get_instrument_data_to_process();
    $self->status_message(sprintf('Found %s instrument data to process.', scalar(@instrument_data)));

    for my $current_inst_data (@instrument_data) {
        my $analysis_project = $self->_get_analysis_project_for_instrument_data($current_inst_data);
        my $config = $analysis_project->get_configuration_reader();
        my $hashes = $self->_prepare_configuration_hashes_for_instrument_data($current_inst_data, $config);
        for my $model_type (keys %$hashes) {
            for my $model_instance (@{$hashes->{$model_type}}) {
                my ($model, $created_new) = $self->_get_model_for_config_hash($model_type, $model_instance);
                $self->status_message(sprintf('Model: %s %s for instrument data: %s.',
                        $model->id, ($created_new ? ' created' : ' found'), $current_inst_data->id ));
                my $success = $self->_assign_instrument_data_to_model($model, $current_inst_data, $created_new);
                $self->_assign_model_to_analysis_project($analysis_project, $model);
                $self->_mark_instrument_data_status($current_inst_data, $success);
            }
        }
    }

    return 1;
}

sub _assign_instrument_data_to_model {
    my ($self, $model, $instrument_data, $newly_created) = @_;

    #if a model is newly created, we want to assign all applicable instrument data to it
    my %params_hash = (model => $model);
    if ($newly_created) {
        $params_hash{all} = 1;
    } else {
        $params_hash{instrument_data} = [$instrument_data];
    }

    my $cmd = Genome::Model::Command::InstrumentData::Assign->create(%params_hash);
    unless ($cmd->execute()) {
        $self->error_message('Failed to assign ' . $instrument_data->__display_name__ . ' to ' . $model->__display_name__);
        return 0;
    }

    $model->build_requested(1);
    return 1;
}

sub _mark_instrument_data_status {
    my ($self, $instrument_data, $success) = @_;

    $instrument_data->remove_attribute(attribute_label => 'tgi_lims_status');
    $instrument_data->remove_attribute(attribute_label => 'tgi_lims_fail_message');

    if ($success) {
        $self->_mark_instrument_data_as_processed($instrument_data);
    } else {
        $self->_mark_instrument_data_as_failed($instrument_data);
    }

    return 1;
}

sub _mark_instrument_data_as_processed {
    my ($self, $instrument_data) = @_;

    $instrument_data->add_attribute(
        attribute_label => 'tgi_lims_status',
        attribute_value => 'processed',
    );
    $instrument_data->remove_attribute(attribute_label => 'tgi_lims_fail_count');

    return 1;
}

sub _mark_instrument_data_as_failed {
    my ($self, $instrument_data) = @_;

    my $fail_count_attr = $instrument_data->attributes(attribute_label => 'tgi_lims_fail_count');
    my $previous_count = 0;
    if ($fail_count_attr) {
        $previous_count = $fail_count_attr->attribute_value;
        $fail_count_attr->delete;
    }

    $instrument_data->add_attribute(
        attribute_label => 'tgi_lims_fail_count',
        attribute_value => ($previous_count + 1),
    );

    $instrument_data->add_attribute(
        attribute_label => 'tgi_lims_status',
        attribute_value => 'failed',
    );

    $instrument_data->add_attribute(
        attribute_label => 'tgi_lims_fail_message',
        attribute_value => $self->error_message,
    );

    return 1;
}

sub _get_model_for_config_hash {
    my $self = shift;
    my $class_name = shift;
    my $config = shift;

    my $m = $class_name->get(%$config);
    #return the model, plus a 'boolean' value indicating if we created a new model
    my @model_info =  $m ? ($m, 0) : ($class_name->create(%$config), 1);
    return wantarray ? @model_info : $model_info[0];
}

sub _prepare_configuration_hashes_for_instrument_data {
    my ($self, $inst_data, $config_obj) = @_;
    #TODO eventually this will need to support multiple references
    my $config_hash = $config_obj->get_config(
        sequencing_platform => $inst_data->sequencing_platform,
        domain              => _normalize_domain($inst_data->taxon->domain),
        taxon               => _normalize_taxon($inst_data->species_name),
        type                => _normalize_extraction_type($inst_data->sample->extraction_type),
    );
    for my $model_type (keys %$config_hash) {
        if (ref $config_hash->{$model_type} ne 'ARRAY') {
            $config_hash->{$model_type} = [$config_hash->{$model_type}];
        }

        for my $model_instance (@{$config_hash->{$model_type}}) {
            my $instrument_data_properties = delete $model_instance->{instrument_data_properties};
            if($instrument_data_properties) {
                while((my $model_property, my $instrument_data_property) = each %$instrument_data_properties) {
                    $model_instance->{$model_property} = $inst_data->$instrument_data_property;
                }
            }
        }
    }
    return $config_hash;
}

sub _get_instrument_data_to_process {
    my $self = shift;

    if ($self->instrument_data) {
        return $self->instrument_data;
    } else {
        return Genome::InstrumentData->get(
            'tgi_lims_status' => [qw/ new failed /],
            'analysis_project_id !=' => undef,
            -hint => [ 'sample', 'sample.source', 'sample.source.taxon', ],
        );
    }
}

sub _get_analysis_project_for_instrument_data {
    my $self = shift;
    my $instrument_data = shift;
    die('You must provide a single piece of instrument data!') unless $instrument_data;

    my $analysis_project_id = $instrument_data->analysis_project_id;
    my $analysis_project = Genome::Config::AnalysisProject->get($analysis_project_id)
        || die("$analysis_project_id doesn't seem to be a valid Genone::Config::AnalysisProject ID!");
    return $analysis_project;
}

sub _assign_model_to_analysis_project {
    my $self = shift;
    my $analysis_project = shift;
    my $model = shift;

    die('Must specify an anlysis project and a model!') unless $analysis_project && $model;

    return $analysis_project->model_group->assign_models($model);
}

sub _lock {
    my $lock_var = $ENV{GENOME_LOCK_DIR} . '/genome_config_command_configure-queued-instrument-data/lock';
    my $lock = Genome::Sys->lock_resource(resource_lock => $lock_var, max_try => 1);

    die('Unable to acquire the lock! Is ConfigureQueuedInstrumentData already running or did it exit uncleanly?')
        unless $lock;

    UR::Context->current->add_observer(
        aspect => 'commit',
        callback => sub {
            Genome::Sys->unlock_resource(resource_lock => $lock);
        }
    );
}

#It's lame that these methods need to exist - is there a way to clean the data before it comes over?
sub _normalize_domain {
    my $domain = lc(shift);
    if ($domain eq 'unknown') {
        return 'eukaryota';
    }
    return $domain;
}

sub _normalize_taxon {
    my $taxon = lc(shift);
    if ($taxon =~ /mus musculus/i) {
        return 'mus_musculus';
    } elsif ($taxon =~ /homo sapien/i) {
        return 'homo_sapiens';
    } elsif ($taxon =~ /maize/i) {
        return 'maize';
    }
}

sub _normalize_extraction_type {
    my $type = lc(shift);
    if ($type =~ /rna/i || $type =~ /cdna/i) {
        return 'rna';
    } elsif ($type =~ /[^c]dna/i) {
        return 'genomic_dna';
    }
    return $type
}

1;
