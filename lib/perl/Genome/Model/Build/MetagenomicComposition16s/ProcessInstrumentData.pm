package Genome::Model::Build::MetagenomicComposition16s::ProcessInstrumentData;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::MetagenomicComposition16s::ProcessInstrumentData {
    is => 'Command::V2',
    has_input => [
        build => { is => 'Genome::Model::Build::MetagenomicComposition16s',
            is_output => 1},
    ],
};

sub execute {
    my $self = shift;

    my $process_ok;
    if ( $self->build->sequencing_platform eq 'sanger' ) {
        my $cmd = Genome::Model::Build::MetagenomicComposition16s::ProcessSangerInstrumentData->create(
            build => $self->build,
        );
        $process_ok = $cmd->prepare_instrument_data;
    }
    else {
        $process_ok = $self->build->prepare_instrument_data;
    }

    if ( not $process_ok ) {
        $self->error_message('Failed to prepare instrument data for '.$self->build->description);
        return;
    }

    return 1;
}

1;

