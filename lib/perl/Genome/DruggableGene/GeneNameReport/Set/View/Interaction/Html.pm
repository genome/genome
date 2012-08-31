package Genome::DruggableGene::GeneNameReport::Set::View::Interaction::Html;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::GeneNameReport::Set::View::Interaction::Html {
    is => 'Genome::View::Status::Html',
    has => {
        perspective => { is => 'Text', value => 'interaction' },
    },
    has_optional => [
        data => { is => 'HASH' },
    ],
};

sub _get_xml_view {
    my $self = shift;
    return Genome::DruggableGene::GeneNameReport::Set::View::Interaction::Xml->create(
        data => $self->data,
        subject => $self->subject,
    );
}

1;
