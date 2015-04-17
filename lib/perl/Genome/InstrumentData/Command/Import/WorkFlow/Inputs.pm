package Genome::InstrumentData::Command::Import::WorkFlow::Inputs;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Import::WorkFlow::Inputs { 
    is => 'UR::Object',
    has => {
        _inputs => { is => 'HASH', },
    },
};

sub create {
    my ($class, %params) = @_;

    my %inputs;
    $inputs{instrument_data_properties} = $class->_resolve_instrument_data_properties($params{instrument_data_properties});

    return $class->SUPER::create(_inputs => \%inputs);
}

sub for_worklflow {
    my $self = shift;
    return $self->_inputs;
}

sub instrument_data_properties {
    my $self = shift;
    return $self->_inputs->{instrument_data_properties};
}

sub _resolve_instrument_data_properties {
    my ($class, $incoming_properties) = @_;

    return {} if not $incoming_properties and not @$incoming_properties;

    my %properties;
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

