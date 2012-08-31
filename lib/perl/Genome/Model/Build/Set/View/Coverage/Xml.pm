package Genome::Model::Build::Set::View::Coverage::Xml;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Set::View::Coverage::Xml {
    is => 'Genome::Model::Set::View::Coverage::Xml',
};

sub members {
    my $self = shift;

    my @b = $self->builds;

    return Genome::Model->get([map($_->model_id, @b)]);
}

sub builds {
    my $self = shift;

    my $set = $self->subject;
    return $set->members;
}

1;
