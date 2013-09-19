package Genome::Workflow::Detail::Operation;

use strict;
use warnings;

use Genome;
use Genome::Workflow::Detail::TypeMap;
use IO::Scalar qw();
use Set::Scalar qw();
use XML::LibXML qw();
use Carp qw(confess);


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

sub from_xml_element {
    my ($class, $element) = @_;

    # Prevent accidental recursion when subclasses don't override this method
    unless ($class eq 'Genome::Workflow::Detail::Operation') {
        confess $class->error_message(sprintf(
                "from_xml_element not implemented in subclass %s", $class));
    }

    my $subclass = $class->_get_subclass_from_element($element);
    return $subclass->from_xml_element($element);
}

sub input_properties {}
sub output_properties {}
sub operation_type_attributes {}


# ------------------------------------------------------------------------------
# Public methods
# ------------------------------------------------------------------------------

sub from_xml {
    my ($class, $xml) = @_;
    my $fh = new IO::Scalar \$xml;

    return $class->from_xml_file($fh);
}

sub from_xml_file {
    my ($class, $fh) = @_;
    my $doc = XML::LibXML->load_xml(IO => $fh);
    return $class->from_xml_element($doc->documentElement);
}

sub from_xml_filename {
    my ($class, $filename) = shift;

    my $fh = Genome::Sys->open_file_for_reading($filename);
    return $class->from_xml_file($fh);
}

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

sub _get_subclass_from_element {
    my ($class, $element) = @_;
    my $nodes = $element->find('operationtype');
    my $operation_type_element = $nodes->pop;
    return Genome::Workflow::Detail::TypeMap::class_from_type(
        $operation_type_element->getAttribute('typeClass'));
}

1;
