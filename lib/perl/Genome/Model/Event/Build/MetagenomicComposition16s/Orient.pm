package Genome::Model::Event::Build::MetagenomicComposition16s::Orient;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::MetagenomicComposition16s::Orient {
    is => 'Genome::Model::Event::Build::MetagenomicComposition16s',
};

sub execute {
    my $self = shift;

    unless ( $self->build->orient_amplicons ) {
        $self->error_message("Failed to orient amplicons for ".$self->build->description);
        return;
    }

    return 1;
}

1;

