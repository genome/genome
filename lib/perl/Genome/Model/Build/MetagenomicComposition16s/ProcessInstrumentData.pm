package Genome::Model::Build::MetagenomicComposition16s::ProcessInstrumentData;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::MetagenomicComposition16s::ProcessInstrumentData {
    is => 'Command::V2',
    has_input => [
        build => { 
            is => 'Genome::Model::Build::MetagenomicComposition16s',
            is_output => 1,
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
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

