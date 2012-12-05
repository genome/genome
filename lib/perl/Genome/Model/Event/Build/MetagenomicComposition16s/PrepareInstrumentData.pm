package Genome::Model::Event::Build::MetagenomicComposition16s::PrepareInstrumentData;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::MetagenomicComposition16s::PrepareInstrumentData {
    is => 'Genome::Model::Event::Build::MetagenomicComposition16s',
};

sub bsub {
    return "-R 'span[hosts=1] select[type=LINUX64]'";
}

sub execute {
    my $self = shift;

    my $process_ok;
    if ( $self->sequencing_platform eq 'sanger' ) {
        my $cmd = Genome::Model::MetagenomicComposition16s::Command::ProcessSangerInstrumentData->create(
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

