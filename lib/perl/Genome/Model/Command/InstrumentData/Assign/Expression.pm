package Genome::Model::Command::InstrumentData::Assign::Expression;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::InstrumentData::Assign::Expression {
    is => 'Genome::Model::Command::InstrumentData::Assign::Base',
    has_input => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            doc => 'instrument data to assign (resolved by expression)',
            shell_args_position => 1,
        },
    ],
    has_optional_input => [
        filter => {
            is => 'Text',
            valid_values => ['forward-only','reverse-only'],
        },
    ],
};

sub help_brief {
    return "Assign specific instrument data to a model by boolean expression";
}

sub help_detail {
    return <<'EOHELP'
This command assigns specific instrument data to the model.

The specification can be a list of IDs or a UR boolean expression.
e.g.:
123456,123457,123458
library.name=H_EX-ample-123_lib456
EOHELP
;
}

sub _resolve_instrument_data {
    my $self = shift;

    return $self->instrument_data;
}

sub _resolve_filter_for_instrument_data {
    my $self = shift;
    my $instrument_data = shift;

    return $self->filter;
}

1;
