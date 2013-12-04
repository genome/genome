package Genome::Config::AnalysisProject;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject {
    is => ['Genome::Utility::ObjectWithTimestamps', 'Genome::Utility::ObjectWithCreatedBy'],
    id_generator => '-uuid',
    data_source => 'Genome::DataSource::GMSchema',
    table_name => 'config.analysis_project',
    id_by => [
        id => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        _model_group_id => {
            is => 'Text',
            column_name => 'model_group_id',
        },
        model_group => {
            is => 'Genome::ModelGroup',
            id_by => '_model_group_id',
        },
        name => {
            is => 'Text',
        },
        status => {
            is => 'Text',
            default_value => 'Pending',
            valid_values => ['Pending', 'Approved', 'In Progress', 'Completed', 'Archived'],
        },
        subject_pairings => {
            is_many => 1,
            is => 'Genome::Config::AnalysisProject::SubjectPairing',
            reverse_as => 'analysis_project',
        },
        analysis_project_bridges => {
            is => 'Genome::Config::AnalysisProject::InstrumentDataBridge',
            reverse_as => 'analysis_project',
            is_many => 1,
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            via => 'analysis_project_bridges',
            to => 'instrument_data',
            is_many => 1,
        },
        samples => {
            is => 'Genome::Subject',
            via => 'instrument_data',
            to => 'sample',
        },
        models => {
            is => 'Genome::Model',
            to => 'models',
            via => 'model_group',
            is_many => 1,
        },
        config_items => {
            is => 'Genome::Config::Profile::Item',
            is_many => 1,
            reverse_as => 'analysis_project',
        },
    ],
    has_transient_optional => [
        configuration_profile => {
            is => 'Genome::Config::Profile'
        }
    ],
};

sub __display_name__ {
    my $self = shift;
    return sprintf('%s (%s)', $self->name, $self->id);
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    eval {
        $self->_create_model_group();
    };
    if(my $error = $@) {
        $self->delete();
        die($error);
    }
    return $self;
}

sub delete {
    my $self = shift;
    eval {
        if ($self->model_group) {
            $self->model_group->delete();
        }
    };
    if(my $error = $@) {
        die($error);
    }
    return $self->SUPER::delete();
}

sub get_configuration_profile {
    my $self = shift;

    unless ($self->configuration_profile) {
        $self->configuration_profile(
            Genome::Config::Profile->create_from_analysis_project($self)
        );
    }

    return $self->configuration_profile;
}

sub _create_model_group {
    my $self = shift;
    my $mg_name = sprintf("%s - %s - Analysis Project", $self->name, $self->id);
    my $mg = Genome::ModelGroup->create(name => $mg_name);
    $self->model_group($mg);
}

1;
