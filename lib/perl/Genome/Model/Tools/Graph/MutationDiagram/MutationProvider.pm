package Genome::Model::Tools::Graph::MutationDiagram::MutationProvider;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Graph::MutationDiagram::MutationProvider {
};

sub next {
    my $self = shift;
    Carp::confess $self->error_message("Must define subroutine next in subclass");
}

1;

