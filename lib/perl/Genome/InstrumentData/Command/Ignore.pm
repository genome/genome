package Genome::InstrumentData::Command::Ignore;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Ignore {
    is => ['Command::V2'],
    has_input => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'The InstrumentData to be modified.',
            is_many => 1,
        },
        ignore => {
            is => 'Number',
            doc => 'Whether to ignore or unignore the given InstrumentData. 1 for ignore, 0 for unignore.'
        },
    ],
};

sub execute {
    my $self = shift;
    my $ignore = $self->ignore();
    if($ignore != 1 and $ignore != 0) {
        $self->error_message("Ignore is not 1 or 0. (Got: $ignore).");
        return 0;
    }
    foreach my $instrument_data ($self->instrument_data) {
        $self->status_message("Trying instrument-data id: " . $instrument_data->id);
        $instrument_data->ignored($self->ignore());
    }
    return 1;
}
