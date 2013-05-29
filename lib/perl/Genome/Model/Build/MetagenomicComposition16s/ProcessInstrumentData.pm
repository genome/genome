package Genome::Model::Build::MetagenomicComposition16s::ProcessInstrumentData;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::MetagenomicComposition16s::ProcessInstrumentData {
    is => 'Command::V2',
    has_input => [
        input_build => {
            is => 'Genome::Model::Build::MetagenomicComposition16s',
            is_many => 1,
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
        },
    ],
    has_output => [
        build => {
            is => 'Genome::Model::Build::MetagenomicComposition16s',
            calculate_from => ['input_build'],
            calculate => sub { return $_[0]; }

        },
    ],
};

sub execute {
    my $self = shift;

    my $process_ok = $self->build->process_instrument_data( $self->instrument_data );
    if ( not $process_ok ) {
        $self->error_message('Failed to process instrument data for '.$self->build->description);
        return;
    }

    return 1;
}

1;

