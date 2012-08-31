package Genome::Model::Build::View::Coverage::Xml;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::View::Coverage::Xml {
    is => 'Genome::Model::Build::Set::View::Coverage::Xml',
};

sub builds {
    my $self = shift;

    return $self->subject;
}

1;
