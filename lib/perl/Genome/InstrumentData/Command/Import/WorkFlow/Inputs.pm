package Genome::InstrumentData::Command::Import::WorkFlow::Inputs;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Import::WorkFlow::Inputs { 
    is => 'UR::Object',
    has_transient => {
        incoming_params => { is => 'HASH', },
        instrument_data_properties => { is => 'HASH', },
        source_files => { is =>'ARRAY', },
    },
};

sub create {
    my ($class, %params) = @_;

    my $self = $class->SUPER::create(incoming_params => \%params);
    return if not $self;

    $self->_resolve_source_files;
    $self->_resolve_instrument_data_properties;

    return $self;
}

sub for_worklflow {
    my $self = shift;
    return {
        source_files => $self->source_files,
        instrument_data_properties => $self->instrument_data_properties,
    };
}

sub _resolve_source_files {
    my $self = shift;

    my $source_files = delete $self->incoming_params->{source_files};
    die $self->error_message('No source files!') if not $source_files;
    my @source_files = split(/,/, $source_files);
    # FIXME check exists?
    
    return $self->source_files(\@source_files);
}

sub _resolve_instrument_data_properties {
    my $self = shift;

    my $incoming_params = $self->incoming_params;
    my $incoming_properties = delete $incoming_params->{instrument_data_properties} || [];
    for my $name (qw/ description downsample_ratio /) {
        my $value = delete $incoming_params->{$name};
        next if not $value;
        push @$incoming_properties, $name.'='.$value;
    }

    my $properties = $self->_resolve_incoming_instrument_data_property_strings($incoming_properties);

    if ( not $properties->{original_data_path} ) {
        $properties->{original_data_path} = join(',', @{$self->source_files});
    }

    return $self->instrument_data_properties($properties);
}

sub _resolve_incoming_instrument_data_property_strings {
    my ($class, $incoming_properties) = @_;

    my %properties;
    return \%properties if not $incoming_properties and not @$incoming_properties;

    for my $key_value_pair ( @$incoming_properties ) {
        my ($label, $value) = split('=', $key_value_pair);
        if ( not defined $value or $value eq '' ) {
            die $class->error_message('Failed to parse with instrument data property label/value! '.$key_value_pair);
        }
        if ( exists $properties{$label} and $value ne $properties{$label} ) {
            die $class->error_message(
                "Multiple values for instrument data property! $label => ".join(', ', sort $value, $properties{$label})
            );
        }
        $properties{$label} = $value;
    
    }

    return \%properties;
}

1;

