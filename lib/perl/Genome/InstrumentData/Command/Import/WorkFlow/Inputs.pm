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
    $inputs{source_files} = $class->_resolve_source_files(\%params);
    $inputs{instrument_data_properties} = $class->_resolve_instrument_data_properties(\%params);

    if ( not $inputs{instrument_data_properties}->{original_data_path} ) {
        $inputs{instrument_data_properties}->{original_data_path} = join(',', @{$inputs{source_files}});
    }

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

sub _resolve_source_files {
    my ($self, $params) = @_;

    my $source_files = delete $params->{source_files};
    die $self->error_message('No source files!') if not $source_files;
    my @source_files = split(/,/, $source_files);
    # FIXME check exists?

    return \@source_files;
}

sub _resolve_instrument_data_properties {
    my ($class, $params) = @_;

    my $incoming_properties = delete $params->{instrument_data_properties} || [];
    for my $name (qw/ description downsample_ratio /) {
        my $value = delete $params->{$name};
        next if not $value;
        push @$incoming_properties, $name.'='.$value;
    }

    return $class->_resolve_incoming_instrument_data_property_strings($incoming_properties);
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

