package Genome::InstrumentData::Command::Ignore;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Ignore {
    is => ['Genome::Command::Base'],
    has_input => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'The InstrumentData to be modified.',
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
    if($ignore != 1 and $ignore != 0)
    {
        $self->error_message("Ignore is not 1 or 0. (Got: $ignore).");
        return 0;
    }
    $self->instrument_data->ignored($self->ignore() );
    return 1;
}
