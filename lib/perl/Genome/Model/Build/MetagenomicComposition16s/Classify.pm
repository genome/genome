package Genome::Model::Build::MetagenomicComposition16s::Classify;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::MetagenomicComposition16s::Classify {
    is => 'Command::V2',
    has_input => [
        input_build => {
            is => 'Genome::Model::Build::MetagenomicComposition16s',
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

    unless ( $self->build->classify_amplicons ) {
        $self->error_message("Failed to classify amplicons for ".$self->build->description);
        return;
    }

    return 1;
}

1;

