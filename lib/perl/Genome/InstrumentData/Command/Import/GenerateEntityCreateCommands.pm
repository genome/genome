package Genome::InstrumentData::Command::Import::GenerateEntityCreateCommands;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Command::Import::CsvParser;
require File::Basename;
require List::MoreUtils;
use Params::Validate qw( :types );

class Genome::InstrumentData::Command::Import::GenerateEntityCreateCommands {
    is => 'Command::V2',
    has_input => {
        file => {
            is => 'Text',
            doc => Genome::InstrumentData::Command::Import::CsvParser->csv_help,
        },
    },
    has_optional_output => {
        output_file => {
            is => 'Text',
            default_value => '-',
            doc => 'Output file to put the commands to create the needed entities for import.',
        },
    },
    has_optional_transient => {
        _names_seen => { is => 'HASH', default_value => {}, },
        _output => { is => 'ARRAY', default_value => [], },
    },
    doc => 'generate commands to create individuals, samples, and libraries that are necessary to import instrument data',
};

sub help_detail {
    return 'Using a modified metadata spreadsheet saved as a comma/tab separarted file, generate the commands necessary to make susupport entities (individaul, sample, library) to import instrument data.';
}

sub execute {
    my $self = shift;
    
    my $parser = Genome::InstrumentData::Command::Import::CsvParser->create(file => $self->file);
    while ( my $entity_params = $parser->next ) {
        # Add entity commands as needed
        my $library = Genome::Library->get(name => $entity_params->{library}->{name});
        next if $library;

        if ( not Genome::Individual->get(name => $entity_params->{individual}->{name}) ) {
            $self->_add_individual_create_command($entity_params);
        }

        if ( not Genome::Sample->get(name => $entity_params->{sample}->{name}) ) {
            $self->_add_sample_create_command($entity_params);
        }

        $self->_add_library_create_command($entity_params);
    }

    $self->_generate_output;

    return 1;
}

sub _add_individual_create_command {
    my ($self, $entity_params) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => HASHREF});
    my $params = $entity_params->{individual};
    return $self->_add_create_command_for_type('individual', $params);
}

sub _add_sample_create_command {
    my ($self, $entity_params) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => HASHREF});
    my $params = $entity_params->{sample};
    die $self->error_message('Sample extraction_type is required to create a sample!') if not $params->{extraction_type}; # FIXME
    $params->{source} = $entity_params->{individual}->{name}; # add the individual name to the params
    return $self->_add_create_command_for_type('sample', $params);
}

sub _add_library_create_command {
    my ($self, $entity_params) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => HASHREF});
    my $params = $entity_params->{library};
    $params->{sample} = $entity_params->{sample}->{name}; # add the sample name to the params
    return $self->_add_create_command_for_type('library', $params);
}

sub _add_create_command_for_type {
    my ($self, $type, $params, $dependent_params) = Params::Validate::validate_pos(
        @_, {type => HASHREF}, {type => SCALAR}, {type => HASHREF},
    );

    return 1 if $self->_names_seen->{$type}->{ $params->{name} };

    # Load create class
    my $create_class = 'Genome::'.ucfirst($type).'::Command::Create';
    my $meta = $create_class->__meta__;
    my %properties = map { $_->property_name => $_ } $meta->properties;
    my %required_properties = map { $_->property_name => 1 } grep { not $_->is_optional } values %properties;
    delete $required_properties{id};

    # Build command
    my $cmd = "genome $type create";
    for my $property_name ( keys %$params ) {
        my $property = $meta->property_meta_for_name($property_name);
        die $self->error_message('Unknown %s property: %s', $type, $property_name) if not $property;
        my $cli_property_name = $property_name;
        $cli_property_name =~ s/_/-/g;
        $cmd .= ' --'.$cli_property_name.'=';
        if ( $property->data_type =~ /^Genome::/ ) { # assuming name
            $cmd .= '"name='.$params->{$property_name}.'"';
        }
        else {
            $cmd .= '"'.$params->{$property_name}.'"';
        }
        delete $required_properties{$property_name};
    }

    # Check required params
    die $self->error_message('Missing required %s properties: %s', $type, join(' ', keys %required_properties)) if %required_properties;

    $self->_names_seen->{$type}->{ $params->{name} } = 1;
    push @{$self->_output}, $cmd."\n";

    return 1;
}

sub _generate_output {
    my $self = shift;
    my @output = @{$self->_output};
    return 1 if not @output;
    return Genome::Sys->write_file($self->output_file, "#!/usr/bin/env bash\n", @output);
}

1;

