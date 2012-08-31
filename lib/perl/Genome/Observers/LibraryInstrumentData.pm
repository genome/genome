package Genome::Observers::LibraryInstrumentData;

use strict;
use warnings;

use Carp;

Genome::Library->add_observer(
    aspect => 'delete',
    callback => \&delete_instrument_data,
);

sub delete_instrument_data {
    my $self = shift;
    my @instrument_data = Genome::InstrumentData->get(library_id => $self->id);
    if (@instrument_data) {
        Carp::confess "Cannot delete library " . $self->__display_name__ . 
            ", there are " . scalar @instrument_data . " derived from it!";
    }
    return 1;
}

1;

