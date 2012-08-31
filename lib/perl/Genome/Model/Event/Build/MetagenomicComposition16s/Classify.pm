package Genome::Model::Event::Build::MetagenomicComposition16s::Classify;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Event::Build::MetagenomicComposition16s::Classify {
    is => 'Genome::Model::Event::Build::MetagenomicComposition16s',
};

sub execute {
    my $self = shift;

    unless ( $self->build->classify_amplicons ) {
        $self->error_message("Failed to classify amplicons for ".$self->build->description);
        return;
    }

    return 1;
}

1;

