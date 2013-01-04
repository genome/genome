package Genome::Model::Build::MetagenomicComposition16s::MergeProcessedInstrumentData;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::MetagenomicComposition16s::MergeProcessedInstrumentData {
    is => 'Command::V2',
    has_input => [
        build => { 
            is => 'Genome::Model::Build::MetagenomicComposition16s',
            is_output => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    my $merge_ok = $self->build->merge_processed_instrument_data;
    if ( not $merge_ok ) {
        $self->error_message('Failed to merge processed instrument data for '.$self->build->description);
        return;
    }

    return 1;
}

1;

