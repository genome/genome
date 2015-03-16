package Genome::WorkflowBuilder::StreamProcess;

use strict;
use warnings;
use Genome;

class Genome::WorkflowBuilder::StreamProcess {
    has => [
        command_name => {
            is => 'String',
        },
        args => {
            is => 'ARRAY',
        },
        in_file_link => {
            is => 'Path',
            is_optional => 1,
        },
        out_file_link => {
            is => 'Path',
            is_optional => 1,
        },
        err_file_link => {
            is => 'Path',
            is_optional => 1,
        },
    ],
};

sub get_xml_element {
    my $self = shift;
    my $element = XML::LibXML::Element->new('command');
    $element->setAttribute('command name', $self->command_name);
    $element->addChild($self->_get_args_xml);
    return $element;
}

sub get_link_xml_elements {
    my $self = shift;
    my @elements;
    if (defined $self->in_file_link) {
        my $element = XML::LibXML::Element->new('connect_input_file');
        $element->setAttribute('source', $self->in_file_link);
        $element->setAttribute('target', $self->name);
        $element->setAttribute('target_fd', 'stdin');
        push @elements, $element;
    }
    if (defined $self->out_file_link) {
        my $element = XML::LibXML::Element->new('connect_output_file');
        $element->setAttribute('source_fd', 'stdout');
        $element->setAttribute('source', $self->name);
        $element->setAttribute('target', $self->out_file_link);
        push @elements, $element;
    }
    if (defined $self->err_file_link) {
        my $element = XML::LibXML::Element->new('connect_output_file');
        $element->setAttribute('source_fd', 'stderr');
        $element->setAttribute('source', $self->name);
        $element->setAttribute('target', $self->out_file_link);
        push @elements, $element;
    }
    return @elements;
}

sub _get_args_xml {
    my $self = shift;
    my $element = XML::LibXML::Element->new('args');
    for my $arg ($self->args) {
        $element->appendTextChild('arg', $arg);
    }
    return $element;
}

1;

