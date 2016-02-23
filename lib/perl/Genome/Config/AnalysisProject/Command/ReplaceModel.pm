package Genome::Config::AnalysisProject::Command::ReplaceModel;

use strict;
use warnings;

use List::MoreUtils 'uniq';
use Params::Validate ':types';
use YAML;

use Genome;

class Genome::Config::AnalysisProject::Command::ReplaceModel {
    is => 'Command::V2',
    has_input => {
       model => {
            is => 'Genome::Model',
            shell_args_position => 1,
            doc => 'Model to replace on the analysis project.',
       },
       new_profile_item => {
            is => 'Genome::Config::Profile::Item',
            shell_args_position => 2,
            doc => 'The configuration to use to copy the model.',
       },
   },
   has_optional_transient => {
       analysis_project => { via => 'new_profile_item', to => 'analysis_project', },
       new_model => { is => 'Genome::Model', },
   },
};

sub help_brief {
    return 'replace a model';
}

sub help_detail {
    return <<EOS
This command takes a model and new config profile item and copies the model, replacing the different items in the new profile config as compared to the model's existing config. The model is queued for building. Additionally, the builds for the old model should be abandoned. That way they do not appear in a failed state.

Use this command when a model needs to replaced on a current analysis project. This is usually due to a need to switch a model parameter because the current configuration fails to produce a successful build. If a different config is needed for most or all instrument data, a new config should be added and the instrument data reprocessed.

The new config profile must have been previously added to the analysis project. Additionally, the profile status must be 'inactive'.

See analysis project configuration documentation for more details on creating/finding a new config item.

EOS
}

sub execute {
    my $self = shift;
    $self->status_message('Replace models on analysis project...');

    $self->_verify_new_config;
    $self->status_message('New config: %s', $self->new_profile_item->file_path);
    $self->status_message('Old model:  %s', $self->model->__display_name__);
    $self->_process_model;
    $self->status_message('NEW model:  %s', $self->new_model->__display_name__);
    $self->status_message('Please abandon builds for old model!');
    return 1;
}

sub _verify_new_config {
    my $self = shift;

    if ( $self->new_profile_item->status ne 'inactive') {
        $self->fatal_message(
            "Config profile (%s) status (%s) must be 'inactive'.", $self->new_profile_item->id, $self->new_profile_item->status
        );
    }

    if ( not $self->analysis_project->is_current ) {
        $self->fatal_message(
            "Analysis project (%s) status (%s) is not a 'current' status.", $self->analysis_project->id, $self->analysis_project->status
        );
    }

    return 1;
}

sub _process_model {
    my $self = shift;

    my @found_models;
    my $config_profile = Genome::Config::Profile->create_from_config_profile_item($self->new_profile_item);
    for my $instrument_data ( $self->model->instrument_data ) {
        push @found_models, $config_profile->process_models_for_instrument_data($instrument_data);
    }

    # filter out duplicates and our orginal model
    @found_models = grep { $_ ne $self->model } uniq @found_models;
    if ( not @found_models ) {
        $self->fatal_message('No model found or created!');
    }
    elsif ( @found_models > 1 ) {
        $self->fatal_message('Found/created more than one model! %s', join(' ', map { $_->__display_name__ } @found_models))
    }

    return $self->new_model($found_models[0]);
}

1;

