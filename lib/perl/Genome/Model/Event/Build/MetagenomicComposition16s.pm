package Genome::Model::Event::Build::MetagenomicComposition16s;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::MetagenomicComposition16s {
    is => 'Genome::Model::Event',
    is_abstract => 1,
    has => [
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            via => 'model',
        },
    ],
};

sub bsub_rusage {
    return "-R 'span[hosts=1]'";
}

sub create {
    my $class = shift;
    
    my $self = $class->SUPER::create(@_)
        or return;

    unless ( $self->build ) {
        $self->error_message("No build given to create.");
        return;
    }

    unless ( $self->model_id or $self->model ) {
        $self->model_id( $self->build->model_id );
    }

    return $self;
}

1;

