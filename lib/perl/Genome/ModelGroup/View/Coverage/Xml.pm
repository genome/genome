package Genome::ModelGroup::View::Coverage::Xml;

use strict;
use warnings;

class Genome::ModelGroup::View::Coverage::Xml {
    is => 'Genome::Model::Set::View::Coverage::Xml',
};

sub members {
    my $self = shift;

    my $group = $self->subject;
    my @members = $group->models;

    return @members;
}

1;
