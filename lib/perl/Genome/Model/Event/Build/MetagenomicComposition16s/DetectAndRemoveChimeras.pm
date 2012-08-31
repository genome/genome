package Genome::Model::Event::Build::MetagenomicComposition16s::DetectAndRemoveChimeras;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::MetagenomicComposition16s::DetectAndRemoveChimeras {
    is => 'Genome::Model::Event::Build::MetagenomicComposition16s',
};

sub execute {
    my $self = shift;

    unless ( $self->build->detect_and_remove_chimeras ) {
        $self->error_message("Failed to detect and remove chimeras for ".$self->build->description);
        return;
    }

    return 1;
}

1;
