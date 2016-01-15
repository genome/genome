package Genome::Config::AnalysisProject::Command::ReplaceModel;

use strict;
use warnings;

use Params::Validate ':types';
use YAML;

use Genome;

class Genome::Config::AnalysisProject::Command::ReplaceModel {
    is => 'Command::V2',
    has_input => {
       model => {
            is => 'Genome::Model',
            doc => 'Model to replace on the analysis project.',
       },
       new_profile_item => {
            is => 'Genome::Config::Profile::Item',
            doc => 'The configuration to use to copy the model.',
       },
   },
   has_transient => {
       analysis_project => { via => 'new_profile_item', to => 'analysis_project', },
       new_model => { is => 'Genome::Model', },
   },
};

sub help_brief {
    return 'replace a model';
}

sub help_detail {
    return <<EOS
This command takes a model and new config profile item and copies the model, replacing the different items in the new profile config as compared to the model's existing config. The model is queued for building. Additionally, the builds forthe old model shold be abandoned. That way they do not appear in a failed state.

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
    my $overrides = $self->_overrides_for_model;
    $self->_copy_model($overrides);

    $self->status_message('NEW model:  %s', $self->new_model->__display_name__);
    return 1;
}

sub _verify_new_config {
    my $self = shift;

    if ( $self->new_profile_item->status ne 'inactive') {
        die $self->error_message(
            "Config profile (%s) status (%s) must be 'inactive'.", $self->new_profile_item->id, $self->new_profile_item->status
        );
    }

    if ( not $self->analysis_project->is_current ) {
        die $self->error_message(
            "Analysis project (%s) status (%s) is not a 'current' status.", $self->analysis_project->id, $self->analysis_project->status
        );
    }

    return 1;
}

sub _overrides_for_model {
    my $self = shift;

    my $config_for_model = $self->_load_model_config_from_file( $self->model->config_profile_item->file_path );
    my $config_for_new_model = $self->_load_model_config_from_file( $self->new_profile_item->file_path );
    my %overrides;
    for my $key ( keys %$config_for_model ) {
        next if ref $config_for_model->{$key};
        next if defined $config_for_new_model->{$key} and $config_for_model->{$key} == $config_for_new_model->{$key};
        my $key_no_id = $key;
        $key_no_id =~ s/_id//;
        $overrides{$key_no_id} = $config_for_new_model->{$key};
    }

    if ( not %overrides ) {
        die $self->error_message("No overrides found for model! %s\nConfig: %s\nNew Config: %s\n",
            $self->model->__display_name__, Data::Dumper::Dumper($config_for_model), Data::Dumper::Dumper($config_for_new_model),
        );
    }

    return \%overrides;
}

sub _load_model_config_from_file {
    my ($self, $file) = Params::Validate::validate_pos(@_, {isa => __PACKAGE__}, {type => SCALAR});

    Genome::Sys->validate_file_for_reading($file);
    my $config = YAML::LoadFile($file);
    die $self->error_message('Failed to load config from file! %s', $file) if not $config;

    my $models_config = $config->{models};
    die $self->error_message('No models key in config: %s', Data::Dumper::Dumper($config)) if not $models_config;

    my $config_for_model;
    for my $key ( keys %{$config->{models}} ) {
        next if not $self->model->class->isa($key);
        $config_for_model = $config->{models}->{$key};
        if ( ref $config->{models}->{$key} eq 'ARRAY' ) { # not currently supporting multiple models for type in config
            die $self->error_message('Model type (%s) config is an array! %s', $key, Data::Dumper::Dumper($config));
        }
        last;
    }

    return $config_for_model;
}

sub _copy_model {
    my ($self, $overrides) = Params::Validate::validate_pos(@_, {isa => __PACKAGE__}, {type => HASHREF});

    my $copy = Genome::Model::Command::Copy->execute(
        model => $self->model,
        overrides => [ map { join('=', $_, $overrides->{$_}) } keys %$overrides ],
    );
    die $self->error_message("Failed to copy model: %s", $self->model->__display_name__) if not $copy->result;
    my $new_model = $copy->_new_model;
    die $self->error_message("Failed to get new model from copy command!") if not $new_model;
    $self->analysis_project->add_model_bridge(model => $new_model, config_profile_item => $self->new_profile_item);
    $new_model->build_requested(1);
    $self->new_model($new_model);

    return 1;
}

1;

