package Genome::Config::Profile;

use strict;
use warnings;

use List::MoreUtils qw/ any uniq /;

use Genome;
use UR::Util;

use feature qw(switch);

class Genome::Config::Profile {
    is => 'UR::Object',
    is_transactional => 0,
    has => [
        config_rule_maps => {
            is_many => 1,
            is => 'Genome::Config::RuleModelMap'
        },
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
        },
    ]
};


sub create_from_analysis_project {
    my $class = shift;
    my $analysis_project = shift;

    my @config_rule_maps = map {
        Genome::Config::Translator->get_rule_model_map_from_config($_)
    } grep{ $_->status eq 'active' } $analysis_project->config_items;

    return $class->create(
        config_rule_maps => \@config_rule_maps,
        analysis_project => $analysis_project,
    );
}

sub create_from_config_profile_item {
    my ($class, $config_profile_item) = @_;

    $class->fatal_message('Cannot use a non current config item to create profile!') if not $config_profile_item->is_current;

    return $class->create(
        config_rule_maps => [ Genome::Config::Translator->get_rule_model_map_from_config($config_profile_item) ],
        analysis_project => $config_profile_item->analysis_project,
    );
}

sub get_config {
    my $self = shift;
    my $inst_data = shift;

    my @model_hashes = map { $_->models }
        grep { $_->match_and_concretize($inst_data) }
        $self->config_rule_maps;

    return $self->_merge_extra_parameters($self->_merge_model_hashes(@model_hashes));
}

sub _merge_model_hashes {
    my $self = shift;
    my @source_hashes = @_;

    my $destination_hash = {};

    for my $source (@source_hashes) {
        for my $key (keys %$source) {
            $destination_hash->{$key} ||= [];
            push @{$destination_hash->{$key}}, @{$source->{$key}};
        }
    }

    return UR::Util::deep_copy($destination_hash);
}

sub _merge_extra_parameters {
    my $self = shift;
    my $config_hash = shift;

    for my $model_config (values %$config_hash) {
        for (@$model_config) {
             $self->_add_user_if_not_present($_);
        }
    }

    return $config_hash;
}

sub _add_user_if_not_present {
    my $self = shift;
    my $hashref = shift;

    $hashref->{run_as} ||= $self->analysis_project->run_as;

    return 1;
}

sub process_models_for_instrument_data {
    my ($self, $instrument_data) = @_;

    my $hashes = $self->_prepare_configuration_hashes_for_instrument_data($instrument_data);
    my @models;
    for my $model_type (keys %$hashes) {
        my $model_hashes = $hashes->{$model_type};
        push @models, $self->_process_models($instrument_data, $model_type, $model_hashes);
    }

    return uniq @models;
}

sub _prepare_configuration_hashes_for_instrument_data {
    my ($self, $instrument_data) = @_;

    my $config_hash = $self->get_config($instrument_data);
    for my $model_type (keys %$config_hash) {
        if (ref $config_hash->{$model_type} ne 'ARRAY') {
            $config_hash->{$model_type} = [$config_hash->{$model_type}];
        }

        my $subject_mapping_attempts = 0;
        my $subject_mapping_successes = 0;

        my @processed_model_instances;

        MODEL_INSTANCE: for my $model_instance (@{$config_hash->{$model_type}}) {
            my $instrument_data_properties = delete $model_instance->{instrument_data_properties};
            if($instrument_data_properties) {
                if(my $input_data = delete $instrument_data_properties->{input_data}) {
                    while ((my $input_name, my $instrument_data_property) = each %$input_data) {
                        $model_instance->{input_data}{$input_name} = $self->_value_for_instrument_data_property($instrument_data, $instrument_data_property);
                    }
                }

                while((my $model_property, my $instrument_data_property) = each %$instrument_data_properties) {
                    $model_instance->{$model_property} = $self->_value_for_instrument_data_property($instrument_data, $instrument_data_property);
                }
            }

            my $requires_subject_mapping = delete $model_instance->{input_data_requires_subject_mapping};
            if ($requires_subject_mapping) {
                $subject_mapping_attempts++;
                my (@processed) = @{ $self->_process_mapped_samples($instrument_data, [{config_profile_item => $model_instance->{config_profile_item} }]) };
                unless (@processed) {
                    #tags didn't match, etc.
                    $self->debug_message('Failed to map subject into configuration.');
                    next MODEL_INSTANCE;
                }

                if (@processed > 1) {
                    $self->fatal_message('Sorry, we do not currently support making multiple models for config-declared subject mappings.');
                }

                my $processed = $processed[0];
                $subject_mapping_successes++;

                delete $processed->{config_profile_item};

                while (my ($key, $value) = each %$processed) {
                    $model_instance->{input_data}{$key} = $value;
                }
            }

            push @processed_model_instances, $model_instance;
        }

        if ($subject_mapping_attempts and not $subject_mapping_successes) {
            $self->fatal_message('Failed to map subject into any configurations for model type %s.', $model_type);
        }

        $config_hash->{$model_type} = \@processed_model_instances;
        $config_hash->{$model_type} = $self->_process_mapped_samples($instrument_data, $config_hash->{$model_type}) if $model_type->requires_subject_mapping;
    }
    return $config_hash;
}

sub _value_for_instrument_data_property {
    my ($self, $instrument_data, $instrument_data_property) = @_;

    if (ref $instrument_data_property eq 'ARRAY') {
        return [
            grep { defined($_) }
            map { $instrument_data->$_ }
            @$instrument_data_property
        ];
    } else {
        return $instrument_data->$instrument_data_property;
    }
}

sub _process_mapped_samples {
    my ($self, $instrument_data, $model_hashes) = @_;
    die('Must provide an analysis project, a piece of instrument data and a config hash!')
        unless($instrument_data && $model_hashes);

    my @subject_mappings = Genome::Config::AnalysisProject::SubjectMapping->get(
        analysis_project => $self->analysis_project,
        subjects => $instrument_data->sample
    );

    unless (@subject_mappings) {
        die(sprintf('Found no mapping information for %s in project %s for a model type that requires mapping!',
            $instrument_data->__display_name__,
            $self->analysis_project->__display_name__));
    }

    return [ map {
        my $mapping = $_;
        my %tags = map { $_->id => 1 } $mapping->tags;
        map { {
          (map { $_->label => $_->subject } $mapping->subject_bridges),
          (map { $_->key => $_->value } $mapping->inputs),
          %$_
        } } grep { $self->_model_hash_matches_tags($_, \%tags) } @$model_hashes
    } @subject_mappings ];
}

sub _model_hash_matches_tags {
    my ($self, $model_hash, $tag_hash) = @_;

    my @model_hash_tags = $model_hash->{config_profile_item}->tags;
    if(keys %$tag_hash) {
        return any { exists $tag_hash->{$_->id} } @model_hash_tags;
    } else {
        return !@model_hash_tags;
    }
}

sub _process_models {
    my $self = shift;
    my $instrument_data = shift;
    my $model_type = shift;
    my $model_list = shift;

    my @models;
    for my $model_instance (@$model_list) {
        my ($model, $created_new, $config_profile_item) = $self->_get_model_for_config_hash($model_type, $model_instance);

        $self->status_message(sprintf('Model: %s %s for instrument data: %s.',
                $model->id, ($created_new ? 'created' : 'found'), $instrument_data->id ));

        $self->_assign_model_to_analysis_project($model, $config_profile_item, $created_new);
        $self->_assign_instrument_data_to_model($model, $instrument_data, $created_new);
        $self->_update_model($model);
        $self->_request_build_if_necessary($model, $created_new);
        push @models, $model;
    }

    return @models;
}

sub _get_model_for_config_hash {
    my $self = shift;
    my $class_name = shift;
    my $config = shift;

    my $config_profile_item = delete $config->{config_profile_item};
    my %read_config = %$config;
    for my $key (keys %read_config) {
        my $value = $read_config{$key};
        if(ref($value) eq 'ARRAY' and scalar(@$value) == 0) {
            $read_config{$key} = undef;
        }
    }

    my @extra_params = (auto_assign_inst_data => 1);
    my @found_models = $class_name->get(@extra_params, %read_config, analysis_project => $self->analysis_project);
    my @m = grep { $_->analysis_project_bridges->profile_item_id eq $config_profile_item->id } @found_models;

    if (scalar(@m) > 1) {
        die(sprintf("Sorry, but multiple identical models were found: %s", join(',', map { $_->id } @m)));
    };

    #return the model, plus a 'boolean' value indicating if we created a new model
    my @model_info;
    if ($m[0]) {
       @model_info = ($m[0], 0, $config_profile_item);
    } else {
       for my $key (keys %$config) {
            delete $config->{$key} unless defined $config->{$key};
       }
       @model_info = ($class_name->create(@extra_params, %$config), 1, $config_profile_item);
    }

    return wantarray ? @model_info : $model_info[0];
}

sub _assign_model_to_analysis_project {
    my $self = shift;
    my $model = shift;
    my $config_profile_item = shift;
    my $created_new = shift;

    die('Must specify an analysis project and a model!') unless $model && $config_profile_item;

    $self->analysis_project->add_model_bridge(model => $model, config_profile_item => $config_profile_item) if $created_new;
    return 1;
}

sub _assign_instrument_data_to_model {
    my ($self, $model, $instrument_data, $newly_created) = @_;

    #if a model is newly created, we want to assign all applicable instrument data to it
    my %params_hash = (model => $model);
    my $cmd = Genome::Model::Command::InstrumentData::Assign::ByExpression->create(
            model => $model,
            instrument_data => [$instrument_data],
            force => 1, #trust the configuration to know what it's doing
        );
    my $executed_ok = eval{ $cmd->execute };

    unless ($executed_ok) {
        die(sprintf('Failed to assign %s to %s', $instrument_data->__display_name__,
                $model->__display_name__));
    }
}

sub _update_model {
    my ($self, $model) = @_;
    if ($model->can("check_for_updates")) {
        unless ($model->check_for_updates) {
            Carp::confess "Could not update model!";
        }
    }
}

sub _request_build_if_necessary {
    my ($self, $model, $newly_created) = @_;

    my $reason = $newly_created? 'created' : 'processed';

    if($model->build_needed) {
        $model->build_requested(1, "CQID $reason model");
    }
}

1;

