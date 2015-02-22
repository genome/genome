package Genome::InstrumentData::Command::Import::WorkFlow::ResolveInstDataProperties;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Import::WorkFlow::ResolveInstDataProperties { 
    is => 'Command::V2',
    has_input => {
        source => {
            is => 'Text',
            doc => 'Source to import.',
        },
    },
    has_optional_input => {
        instrument_data_properties => {
            is => 'Text',
            is_many => 1,
            doc => 'Name and value pairs to add to the instrument data. Separate name and value with an equals (=) and name/value pairs with a comma (,).',
        },
    },
    has_output => {
        resolved_instrument_data_properties => { is => 'Hash', },
    },
};

sub execute {
    my $self = shift;

    my $resolve_instdata_props = $self->_resolve_instrument_data_properties;
    return if not $resolve_instdata_props;

    return 1;
}

sub _resolve_instrument_data_properties {
    my $self = shift;

    my %properties;
    for my $key_value_pair ( $self->instrument_data_properties ) {
        my ($label, $value) = split('=', $key_value_pair);
        if ( not defined $value or $value eq '' ) {
            $self->error_message('Failed to parse with instrument data property label/value! '.$key_value_pair);
            return;
        }
        if ( exists $properties{$label} and $value ne $properties{$label} ) {
            $self->error_message(
                "Multiple values for instrument data property! $label => ".join(', ', sort $value, $properties{$label})
            );
            return;
        }
        $properties{$label} = $value;
    
    }

    $properties{original_data_path} = $self->source;

    return $self->resolved_instrument_data_properties(\%properties);
}

1;

