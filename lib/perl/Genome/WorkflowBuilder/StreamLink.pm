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
    $element->addAttribute('source', sprintf('"%s"', $self->source->name));
    $element->addAttribute('source_fd', sprintf('"%s"', $self->source_fd));
    $element->addAttribute('target', sprintf('"%s"', $self->target->name));
    $element->addAttribute('target_fd', '"stdin"');
    return $element;
}

1;

