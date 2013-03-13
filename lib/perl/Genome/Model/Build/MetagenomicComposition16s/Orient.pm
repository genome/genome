package Genome::Model::Build::MetagenomicComposition16s::Orient;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::MetagenomicComposition16s::Orient {
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

    unless ( $self->build->orient_amplicons ) {
        $self->error_message("Failed to orient amplicons for ".$self->build->description);
        return;
    }

    return 1;
}

1;

