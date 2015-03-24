package Genome::WorkflowBuilder::StreamLink;

use strict;
use warnings;
use Genome;

class Genome::WorkflowBuilder::StreamLink {
    has => [
        source => {
            is => 'Genome::WorkflowBuilder::StreamProcess',
        },
        source_fd => {
            is => 'String',
            valid_values => [qw(stdout stderr)],
        },
        target => {
            is => 'Genome::WorkflowBuilder::StreamProcess',
        }
    ],
};

sub get_xml_element {
    my $self = shift;
    my $element = XML::LibXML::Element->new('connect');
    $element->setAttribute('source', $self->source->name);
    $element->setAttribute('source_fd', $self->source_fd);
    $element->setAttribute('target', $self->target->name);
    $element->setAttribute('target_fd', 'stdin');
    return $element;
}

1;

