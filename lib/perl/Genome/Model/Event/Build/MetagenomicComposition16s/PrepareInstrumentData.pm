package Genome::Model::Event::Build::MetagenomicComposition16s::PrepareInstrumentData;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::MetagenomicComposition16s::PrepareInstrumentData {
    is => 'Genome::Model::Event::Build::MetagenomicComposition16s',
    #is_abstract => 1,
};

sub bsub {
    return "-R 'span[hosts=1] select[type=LINUX64]'";
}

sub execute {
    my $self = shift;

    unless ( $self->build->prepare_instrument_data ) {
        $self->error_message('Failed to prepare instrument data for '.$self->build->description);
        return;
    }

    return 1;
}

1;

#$HeadURL$
#$Id$
