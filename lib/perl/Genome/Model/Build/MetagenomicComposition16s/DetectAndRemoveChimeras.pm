package Genome::Model::Build::MetagenomicComposition16s::DetectAndRemoveChimeras;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::MetagenomicComposition16s::DetectAndRemoveChimeras {
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

    unless ( $self->build->detect_and_remove_chimeras ) {
        $self->error_message("Failed to detect and remove chimeras for ".$self->build->description);
        return;
    }

    return 1;
}

1;

