package Genome::Model::Build::MetagenomicComposition16s::MergeProcessedInstrumentData;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::MetagenomicComposition16s::MergeProcessedInstrumentData {
    is => 'Command::V2',
    has_input => [
        input_build => {
            is => 'Genome::Model::Build::MetagenomicComposition16s',
        },
        dummy_input => {
            is => 'Integer',
            is_optional => 1,
            is_many => 1,
            doc => 'Unused input to facilitate legacy Workflow flow control.',
        },
    ],
    has_output => [
        build => {
            is => 'Genome::Model::Build::MetagenomicComposition16s',
            via => '__self__',
            to => 'input_build',
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

