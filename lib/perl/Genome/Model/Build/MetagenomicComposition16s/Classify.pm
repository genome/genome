package Genome::Model::Build::MetagenomicComposition16s::Classify;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::MetagenomicComposition16s::Classify {
    is => 'Command::V2',
    has_input => [
        input_build => {
            is => 'Genome::Model::Build::MetagenomicComposition16s',
            is_many => 1,
        },
    ],
    has_output => [
        build => {
            is => 'Genome::Model::Build::MetagenomicComposition16s',
            calculate_from => ['input_build'],
            calculate => sub { return $_[0]; }

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

