package Genome::WorkflowBuilder::Block;

use strict;
use warnings;

use Genome;

class Genome::WorkflowBuilder::Block {
    is => 'Genome::WorkflowBuilder::Detail::Operation',
    has => [
        properties => {
            is_many => 1,
            is => 'Text',
        },
    ],
};

sub get_ptero_builder_task {
    require Ptero::Builder::Detail::Workflow::Task;

    my $self = shift;

    $self->validate;

    my %params = (
        name => $self->name,
        methods => [
            $self->_get_ptero_block_method(),
        ],
    );
    if (defined $self->parallel_by) {
        $params{parallel_by} = $self->parallel_by;
    }
    return Ptero::Builder::Detail::Workflow::Task->new(%params);
}

sub _get_ptero_block_method {
    require Ptero::Builder::Detail::Workflow::Converge;

    my $self = shift;
    my ($output_name) = $self->output_properties;
    return Ptero::Builder::Detail::Workflow::Block->new(
        name => 'block',
        parameters => {
            input_names => [$self->input_properties],
            output_name => $output_name,
        },
    );
}

sub input_properties {
    return $_[0]->properties;
}

sub output_properties {
    return $_[0]->properties;
}

sub from_xml_element {
    my ($class, $element) = @_;

    my @properties = map $_->textContent, $element->findnodes('.//property');

    return $class->create(
        name => $element->getAttribute('name'),
        properties => \@properties,
    );
}

sub is_property {
    my ($self, $name) = @_;

    return Set::Scalar->new($self->properties)->contains($name);
}

sub is_input_property {
    my $self = shift;
    return $self->is_property(@_);
}

sub is_output_property {
    my $self = shift;
    return $self->is_property(@_);
}

sub _add_property_xml_elements {
    my ($self, $element) = @_;

    map $self->_add_property_xml_element($element, 'property', $_), $self->properties;

    return;
}


1;
