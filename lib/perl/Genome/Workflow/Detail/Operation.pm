package Genome::Workflow::Detail::Operation;

use strict;
use warnings;

use Genome;
use Genome::Workflow::Detail::TypeMap;
use Set::Scalar;
use XML::LibXML;


class Genome::Workflow::Detail::Operation {
    is => 'Genome::Workflow::Detail::Element',
    is_abstract => 1,

    has => [
        name => {
            is => 'Text',
        },

        log_dir => {
            is => 'Text',
            is_optional => 1,
        },

        parallel_by => {
            is => 'Text',
            is_optional => 1,
        },
    ],
};


# ------------------------------------------------------------------------------
# Abstract methods
# ------------------------------------------------------------------------------

sub input_properties {}
sub output_properties {}
sub operation_type_attributes {}


# ------------------------------------------------------------------------------
# Public methods
# ------------------------------------------------------------------------------

sub is_input_property {
    my ($self, $property_name) = @_;
    return $self->command->__meta__->properties(property_name => $property_name,
        is_input => 1);
}

sub is_output_property {
    my ($self, $property_name) = @_;
    return $self->command->__meta__->properties(property_name => $property_name,
        is_output => 1);
}

sub is_many_property {
    my ($self, $property_name) = @_;
    return $self->command->__meta__->properties(property_name => $property_name,
        is_many => 1);
}

sub operation_type {
    my $self = shift;

    return Genome::Workflow::Detail::TypeMap::type_from_class($self->class);
}


# ------------------------------------------------------------------------------
# Inherited methods
# ------------------------------------------------------------------------------

sub get_xml {
    my $self = shift;

    my $doc = XML::LibXML::Document->new();
    $doc->setDocumentElement($self->get_xml_element);

    return $doc->toString();
}

sub get_xml_element {
    my $self = shift;

    $self->validate;

    my $element = XML::LibXML::Element->new('operation');
    $element->setAttribute('name', $self->name);

    if (defined($self->parallel_by)) {
        $element->setAttribute('parallelBy', $self->parallel_by);
    }

    $element->addChild($self->_get_operation_type_xml_element);

    return $element;
}

my $_INVALID_NAMES = new Set::Scalar('input connector', 'output connector');
sub validate {
    my $self = shift;

    if ($_INVALID_NAMES->contains($self->name)) {
        die $self->error_message(sprintf("Operation name '%s' is not allowed",
            $self->name));
    }

    return;
}


# ------------------------------------------------------------------------------
# Private Methods
# ------------------------------------------------------------------------------

sub _get_operation_type_xml_element {
    my $self = shift;

    my $element = XML::LibXML::Element->new('operationtype');

    $element->setAttribute('typeClass', $self->operation_type);

    map {$self->_add_property_xml_element($element, 'inputproperty', $_)}
        $self->input_properties;
    map {$self->_add_property_xml_element($element, 'outputproperty', $_)}
        $self->output_properties;

    my %attributes = $self->operation_type_attributes;
    for my $attr_name (keys(%attributes)) {
        $element->setAttribute($attr_name, $attributes{$attr_name});
    }

    return $element;
}

sub _add_property_xml_element {
    my ($self, $element, $xml_tag, $text) = @_;

    my $inner_element = XML::LibXML::Element->new($xml_tag);
    $inner_element->appendText($text);
    $element->addChild($inner_element);

    return;
}

1;
