package Genome::Config::AnalysisProject;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject {
    is => 'Genome::Utility::ObjectWithTimestamps',
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
        _configuration_set_id => {
            is => 'Text',
            column_name => 'configuration_set_id',
        },
        _configuration_set => {
            is => 'Genome::Config::Set',
            id_by => '_configuration_set_id',
        },
        created_by => {
            is => 'Text',
        },
        _analysis_menu_item_id => {
            is => 'Text',
            column_name => 'analysis_menu_item_id',
        },
        _analysis_menu_item => {
            is => 'Genome::Config::AnalysisMenuItem',
            id_by => '_analysis_menu_item_id',
        },
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
        created_at => {
            is => 'Timestamp',
        },
        updated_at => {
            is => 'Timestamp',
        },
    ],
    has_transient_optional => [
        configuration_reader => {
            is => 'Genome::Config::MaskedConfigurationReader',
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    eval {
        $self->_create_configuration_set();
        $self->_populate_created_by();
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
        if ($self->_configuration_set) {
            $self->_configuration_set->delete();
        }
        if ($self->model_group) {
            $self->model_group->delete();
        }
    };
    if(my $error = $@) {
        die($error);
    }
    return $self->SUPER::delete();
}

sub get_configuration_reader {
    my $self = shift;

    unless($self->configuration_reader) {
        $self->configuration_reader(Genome::Config::MaskedConfigurationReader->create(
            config_handler              => Genome::Config::Handler::TreeHandler->create(
                                                base_path => $self->_configuration_set->path
                                           ),
            mask_handler                => Genome::Config::Handler::TreeHandler->create(
                                                base_path => $self->_analysis_menu_item->path
                                           ),
            default_handler             => Genome::Config::Handler::TreeHandler->create(
                                                base_path => $ENV{GENOME_ANALYSIS_PROJECT_DEFAULTS}
                                           ),
            configuration_parser        => Genome::Config::Parser::YAML->create(),
            configuration_copy_strategy => Genome::Config::CopyStrategy::TreeCopy->create(),
        ));
    }
    return $self->configuration_reader;
}

sub _populate_created_by {
    my $self = shift;
    unless ($self->created_by) {
        $self->created_by(Genome::Sys->username);
    }
}

sub _create_configuration_set {
    my $self = shift;
    my $set = Genome::Config::Set->create();
    $self->_configuration_set($set);
}

sub _create_model_group {
    my $self = shift;
    my $mg_name = sprintf("%s - %s - Analysis Project", $self->name, $self->id);
    my $mg = Genome::ModelGroup->create(name => $mg_name);
    $self->model_group($mg);
}

1;
