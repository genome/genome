package Genome::Sample::Command::Import::Manager;

use strict;
use warnings;

use Genome;

use Data::Dumper;
use Genome::Sample::Command::Import;
use IO::File;
use Switch;
use YAML;

class Genome::Sample::Command::Import::Manager {
    is => 'Command::V2',
    doc => 'Manage importing a group of samples including importing instrument data and creating and building models.',
    has => [
        working_directory => {
            is => 'Text',
            doc => 'Directory to read and write.',
        },
    ],
    has_optional => [
        status_only => { 
            is => 'Boolean',
            default_value => 0,
            doc => 'Only print status. Do not create samples, models and run imports.', 
        },
    ],
    has_optional_calculated => [
        sample_csv_file => {
            calculate_from => 'working_directory',
            calculate => sub{ my $working_directory = shift; return $working_directory.'/samples.csv'; },
            doc => 'CSV file of samples and attributes. A column called "name" is required. The name should be dash (-) separated values of the nomenclature, indivdual id and sample id.',
        },
        config_file => {
            calculate_from => 'working_directory',
            calculate => sub{ my $working_directory = shift; return $working_directory.'/config.yaml'; },
        },
    ],
    has_optional_transient => [
        config => { is => 'Hash', },
        namespace => { is => 'Text', },
        samples => { is => 'Hash', },
    ],
};

sub execute {
    my $self = shift;

    my $samples = $self->_load_samples;
    return if not $samples;

    my $make_progress = $self->_make_progress;
    return if not $make_progress;

    $self->_status($samples);

    return 1;
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    if ( not -d $self->working_directory ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ working_directory /],
            desc => 'Working directory does not exist or is not aq directory!',
        );
        return @errors;
    }

    my $config_error = $self->_load_config;
    if ( $config_error ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ config_file /],
            desc => $config_error,
        );
        return @errors;
    }

    my $sample_csv_error = $self->_load_sample_csv_file;
    if ( $sample_csv_error ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ sample_csv_file /],
            desc => $sample_csv_error,
        );
        return @errors;
    }

    return;
}

sub _load_config {
    my $self = shift;

    my $config_file = $self->config_file;
    if ( not -s $config_file ) {
        return 'Config file does not exist! '.$config_file;
    }

    my $config = YAML::LoadFile($config_file);
    if ( not $config ) {
        return 'Failed to load config file! '.$config_file;
    }
    $self->config($config);

    my $nomenclature = $config->{sample}->{nomenclature};
    if ( not $nomenclature ) {
        return 'No nomenclature in config file! '.$config_file;
    }

    my $namespace = Genome::Sample::Command::Import->namespace_for_nomenclature($nomenclature);
    if ( not $namespace ) {
        return 'Could not get namespace from miporter. Please ensure there is a config for the '.$nomenclature.' nomenclature.';
    }
    $self->namespace($namespace);

    return;
}

sub _load_sample_csv_file {
    my $self = shift;

    my $sample_csv_file = $self->sample_csv_file;
    if ( not -s $sample_csv_file ) {
        return 'Sample csv file does not exist! '.$sample_csv_file;
    }

    my $sample_csv_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $sample_csv_file,
        separator => ',',
    );
    if ( not $sample_csv_reader ) {
        return 'Failed to open sample csv! '.$sample_csv_file;
    }

    my %headers_not_found = ( name => 1, source_files => 1, );
    for my $header ( @{$sample_csv_reader->headers} ) {
        delete $headers_not_found{$header};
    }

    if ( %headers_not_found ) {
        return 'No '.join(' ', map { '"'.$_.'"' } keys %headers_not_found).' column in sample csv! '.$self->sample_csv_file;
    }

    my @importer_property_names = Genome::Sample::Command::Import->importer_property_names_for_namespace($self->namespace);
    my %samples;
    while ( my $sample = $sample_csv_reader->next ) {
        my $name = delete $sample->{name};
        $samples{$name} = {
            name => $name,
            source_files => [ split(' ', delete $sample->{source_files}) ],
            importer_params => { name => $name, },
        };
        for my $attr ( sort keys %$sample ) {
            my $value = $sample->{$attr};
            next if not defined $value or $value eq '';
            if ( $attr =~ /^(sample|individual)\./ ) { # is sample/individual/inst data indicated?
                push @{$samples{$name}->{importer_params}->{$1.'_attributes'}}, $attr."=\'".$value."\'";
            }
            elsif ( $attr =~ s/^instrument_data\.// ) { # inst data attr
                push @{$samples{$name}->{'instrument_data_attributes'}}, $attr."=\'".$value."\'";
            }
            elsif ( grep { $attr eq $_ } @importer_property_names ) { # is the attr specified in the importer?
                $samples{$name}->{importer_params}->{$attr} = $value;
            }
            else { # assume sample attribute
                push @{$samples{$name}->{importer_params}->{'sample_attributes'}}, $attr."=\'".$value."\'";
            }
        }
    }
    $self->samples(\%samples);

    return;
}

sub _load_samples {
    my $self = shift;

    my $set_job_status_to_samples = $self->_set_job_status_to_samples;
    return if not $set_job_status_to_samples;

    my $samples = $self->samples;
    my %instrument_data = map { $_->sample->name, $_ } Genome::InstrumentData::Imported->get(
        original_data_path => [ map { join(',', @{$_->{source_files}}) } values %$samples ],
        '-hint' => [qw/ attributes /],
    );

    my ($model_class, $model_params) = $self->_resolve_model_params;
    return if not $model_class;

    for my $name ( keys %$samples ) {
        # Inst Data
        my $instrument_data = $instrument_data{$name};
        $samples->{$name}->{instrument_data} = $instrument_data;
        $samples->{$name}->{instrument_data_file} = eval{
            my $attribute = $instrument_data->attributes(attribute_label => 'bam_path');
            if ( not $attribute ) {
                $attribute = $instrument_data->attributes(attribute_label => 'archive_path');
            }
            return $attribute->attribute_value if $attribute;
        };

        # Sample
        my $genome_sample = Genome::Sample->get(name => $name);
        if ( not $genome_sample and not $self->status_only ) {
            $genome_sample = $self->_create_sample($samples->{$name}->{importer_params});
            return if not $genome_sample;
        }
        $samples->{$name}->{sample} = $genome_sample;
        $model_params->{subject} = $genome_sample;

        # Model
        my $model = $model_class->get(%$model_params) if $genome_sample;
        if ( not $model and not $self->status_only ) {
            $model = $model_class->create(%$model_params);
            return if not $model
        }
        $samples->{$name}->{model} = $model;
        if ( $model ) {
            # Model should have instrument data
            if ( $instrument_data ) {
                my @model_instrument_data = $model->instrument_data;
                if ( @model_instrument_data and not grep { $instrument_data->id eq $_ } map { $_->id } @model_instrument_data ) {
                    $model->add_instrument_data($instrument_data);
                }
                # Set latest build
                $samples->{$name}->{build} = $model->latest_build;
            }
        }

        # Import Command
        my $cmd = $self->_resolve_import_command_for_sample($samples->{$name});
        return if not $cmd;
        $samples->{$name}->{import_command} = $cmd;

        # Status
        $self->set_sample_status($samples->{$name});
    }

    $self->samples($samples);
    return 1;
}

sub _set_job_status_to_samples {
    my $self = shift;

    my $samples = $self->samples;
    Carp::confess('Need samples to load job status!') if not $samples;

    return 1 if not $self->config->{'job dispatch'};

    my $job_list_cmd = $self->config->{'job dispatch'}->{list}->{'command'};
    if ( not $job_list_cmd ) {
        $self->error_message('No job list "command" in config! '.YAML::Dump($self->config));
        return;
    }

    my $name_column = $self->config->{'job dispatch'}->{list}->{'name column'};
    if ( not $name_column ) {
        $self->error_message('No job list "name column" in config! '.YAML::Dump($self->config));
        return;
    }
    $name_column--;

    my $status_column = $self->config->{'job dispatch'}->{list}->{'status column'};
    if ( not $status_column ) {
        $self->error_message('No job list "status column" in config! '.YAML::Dump($self->config));
        return;
    }
    $status_column--;

    $job_list_cmd .= ' 2>/dev/null |';
    my $fh = IO::File->new($job_list_cmd);
    while ( my $line = $fh->getline ) {
        my @tokens = split(/\s+/, $line);
        my $name = $tokens[ $name_column ];
        next if not $samples->{$name};
        $samples->{$name}->{job_status} = lc $tokens[ $status_column ];
    }

    return 1;
}

sub _create_sample {
    my ($self, $importer_params) = @_;

    Carp::confess('No params to create sample!') if not $importer_params;

    local $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 0; # quiet the importer
    my $importer_class_name = Genome::Sample::Command::Import->importer_class_name_for_namespace($self->namespace);
    my $importer = $importer_class_name->create(%$importer_params);
    if ( not $importer ) {
        $self->error_message('Failed to create sample importer for sample! '.Data::Dumper::Dumper($importer_params));
        return;
    }
    if ( not $importer->execute ) {
        $self->error_message('Failed to execute sample importer for sample! '.Data::Dumper::Dumper($importer_params));
        return;
    }
    my $sample = $importer->_sample;
    if ( not $sample ) {
        $self->error_message('Executed the importer successfully, but failed to create sample!');
        return;
    }

    return $sample;
}

sub set_sample_status {
    my ($self, $sample) = @_;

    Carp::confess('No sample to set status!') if not $sample;

    $sample->{status} = eval{
        return 'sample_needed' if not $sample->{sample};
        return 'import_'.$sample->{job_status} if $sample->{job_status};
        return 'import_needed' if not $sample->{instrument_data};
        return 'import_failed' if not defined $sample->{instrument_data_file} or not -s $sample->{instrument_data_file};
        return 'model_needed' if $sample->{instrument_data} and not $sample->{model};
        return 'build_requested' if $sample->{model}->build_requested;
        return 'build_needed' if not $sample->{build};
        return 'build_'.lc($sample->{build}->status);
    };

    return $sample->{status};
}

sub _status {
    my $self = shift;

    my $samples = $self->samples;
    Carp::confess('No samples to display status!') if not $samples;

    my %totals;
    my $status;
    for my $name ( sort { $a cmp $b } keys %$samples ) {
        my $sample = $samples->{$name};
        $totals{total}++;
        $self->set_sample_status($sample);
        $totals{ $sample->{status} }++;
        $totals{build}++ if $sample->{status} =~ /^build/;
        $status .= sprintf("%-20s %10s\n", $sample->{name}, $sample->{status});
    }
    print "$status\nSummary:\n".join("\n", map { sprintf('%-16s %s', $_, $totals{$_}) } sort { $a cmp $b } keys %totals)."\n";
    return 1;
}

sub _resolve_model_params {
    my $self = shift;

    my $config = $self->config;
    if ( not $config->{model} ) {
        $self->error_message('Cannot create models. These is no model info in config! '.$config);
        return;
    }

    my %model_params;
    for my $name ( keys %{$config->{model}} ) {
        $model_params{$name} = $config->{model}->{$name};
    }

    if ( not $model_params{processing_profile_id} ) {
        $self->error_message('No processing profile id for model in config! '.Data::Dumper::Dumper($config->{model}));
        return;
    }
    my $processing_profile_id = delete $model_params{processing_profile_id};
    $model_params{processing_profile} = Genome::ProcessingProfile->get($processing_profile_id);
    if ( not $model_params{processing_profile} ) { 
        $self->error_message('Failed to get processing profile for id! '.$processing_profile_id);
        return;
    }

    my $model_class = $model_params{processing_profile}->class;
    $model_class =~ s/ProcessingProfile/Model/;
    my $model_meta = $model_class->__meta__;
    if ( not $model_meta ) {
        $self->error_message('Failed to get model class meta! '.$model_class);
        return;
    }

    for my $param_name ( keys %model_params ) {
        next if $param_name !~ s/_id$//;
        my $param_value_id = delete $model_params{$param_name.'_id'};
        my $param_property = $model_meta->property_meta_for_name($param_name);
        Carp::confess('Failed to get property for model param! '.$param_name) if not $param_property;
        my $param_value_class = $param_property->data_type;
        my $param_value = $param_value_class->get($param_value_id);
        Carp::confess("Failed to get $param_value_class for id! $param_value_id") if not $param_value;
        $model_params{$param_name} = $param_value;
    }

    return ($model_class, \%model_params);
}

sub _resolve_import_command_for_sample {
    my ($self, $sample) = @_;

    Carp::confess('No sample to resolve instrument data import command!') if not $sample;

    my $cmd;
    if ( $self->config->{'job dispatch'}->{launch} ) {
        my $sample_name_cnt = $self->config->{'job dispatch'}->{launch} =~ /\%s/;
        if ( not $sample_name_cnt ) {
            $self->error_message('No "%s" in job dispatch config to replace with sample name.');
            return;
        }
        $cmd = sprintf(
            $self->config->{'job dispatch'}->{launch},
            $sample->{name}, $sample->{name}
        );
        $cmd .= ' ';
    }

    $cmd .= sprintf(
        'genome instrument-data import basic --sample name=%s --source-files %s --import-source-name %s',
        $sample->{name},
        join(',', @{$sample->{source_files}}),
        $self->config->{'sample'}->{'nomenclature'},
    );
    if ( $sample->{instrument_data_attributes} ) {
        $cmd .= ' --instrument-data-properties '.join(',', @{$sample->{instrument_data_attributes}});
    }

    return $cmd;
}

sub _make_progress {
    my $self = shift;

    return 1 if $self->status_only;

    my $samples = $self->samples;
    Carp::confess('No samples to display status!') if not $samples;

    for my $sample ( @$samples ) {
        if ( $sample->{status} eq 'import_needed' ) {
            system $sample->{import_command};
            $self->set_sample_status($sample);
        }
        elsif ( $sample->{status} eq 'import_failed' ) {
            #defined $sample->{instrument_data_file} and -s $sample->{instrument_data_file} );
            # clenaup, relaunch
            #$sample->{inst_data}->delete;
            $self->set_sample_status($sample);
        }
        elsif ( $sample->{status} eq 'build_needed'
            #or $sample->{status} eq 'build_failed'
        ) {
            $sample->{model}->build_requested(1);
            $self->set_sample_status($sample);
        }
    }

    return 1;
}

1;

