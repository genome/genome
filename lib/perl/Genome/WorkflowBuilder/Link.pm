package Genome::WorkflowBuilder::Link;

use strict;
use warnings;

use XML::LibXML qw();


class Genome::WorkflowBuilder::Link {
    is => 'Genome::WorkflowBuilder::Detail::Element',

    has_optional => [
        source => {
            is => 'Genome::WorkflowBuilder::Detail::Operation',
        },

        destination => {
            is => 'Genome::WorkflowBuilder::Detail::Operation',
        },
    ],

    has => [
        source_property => {
            is => 'Text',
        },

        destination_property => {
            is => 'Text',
        },
    ],
};


# ------------------------------------------------------------------------------
# Inherited Methods
# ------------------------------------------------------------------------------

sub get_xml {
    my $self = shift;
    return $self->get_xml_element->toString();
}

sub get_xml_element {
    my $self = shift;

    $self->validate;

    my $element = XML::LibXML::Element->new('link');
    $element->setAttribute('fromOperation', $self->source_operation_name);
    $element->setAttribute('fromProperty', $self->source_property);
    $element->setAttribute('toOperation', $self->destination_operation_name);
    $element->setAttribute('toProperty', $self->destination_property);

    return $element;
}

sub validate {
    my $self = shift;

    $self->_validate_operation_types;

    # make sure x_property is in x as appropriate input/output
    $self->_validate_source_property;
    $self->_validate_output_property;

    # if dest operation is parallel by dest_property, make sure that source_property is_many
    $self->_validate_parallel_by_destination;
}


# ------------------------------------------------------------------------------
# Public Methods
# ------------------------------------------------------------------------------

sub destination_operation_name {
    my $self = shift;
    return $self->_operation_name('destination', 'output connector');
}

sub source_operation_name {
    my $self = shift;
    return $self->_operation_name('source', 'input connector');
}

sub external_input {
    my $self = shift;
    return !defined($self->source);
}

sub external_output {
    my $self = shift;
    return !defined($self->destination);
}

sub sort_key {
    my $self = shift;
    return sprintf("%s|%s|%s|%s",
        $self->_source_name, $self->_destination_name,
        $self->source_property, $self->destination_property);
}


# ------------------------------------------------------------------------------
# Private Methods
# ------------------------------------------------------------------------------

sub _source_name {
    my $self = shift;

    if ($self->source) {
        return $self->source->name;
    } else {
        return 'input connector';
    }
}

sub _destination_name {
    my $self = shift;

    if ($self->destination) {
        return $self->destination->name;
    } else {
        return 'output connector';
    }
}

sub _operation_name {
    my ($self, $operation, $default) = @_;

    if (defined($self->$operation)) {
        return $self->$operation->name;
    } else {
        return $default;
    }
}

sub _validate_operation_types {
    my $self = shift;

    $self->_validate_general_operation_type('source');
    $self->_validate_general_operation_type('destination');

    return;
}

sub _validate_general_operation_type {
    my ($self, $property_name) = @_;

    if (defined($self->$property_name)) {
        unless ($self->$property_name->isa(
                'Genome::WorkflowBuilder::Detail::Operation')) {
            die $self->error_message(sprintf(
                "Expected %s => Genome::WorkflowBuilder::Detail::Operation, "
                . "got '%s' instead", $property_name, $self->$property_name));
        }
    }

}

sub _validate_source_property {
    my $self = shift;

    if (defined($self->source)) {
        unless ($self->source->is_output_property($self->source_property)) {
            die $self->error_message(sprintf(
"Source property '%s' from operation (%s) is not an output",
                    $self->source_property, $self->source->name
            ));
        }
    }
}

sub _validate_output_property {
    my $self = shift;

    if (defined($self->destination)) {
        unless ($self->destination->is_input_property(
                $self->destination_property)) {
            die $self->error_message(sprintf(
"Destination property '%s' from operation (%s) is not an output",
                    $self->destination_property, $self->destination->name
            ));
        }
    }
}

sub _validate_parallel_by_destination {
    my $self = shift;

    if (defined($self->source) && defined($self->destination)) {
        if (defined($self->destination->parallel_by)) {
            unless ($self->source->is_many_property($self->source_property)) {
                die $self->error_message(sprintf(
                    "Source property '%s' (%s) is not is_many for parallel_by "
                    . "destination property '%s' (%s)",
                    $self->source_property, $self->source->name,
                    $self->destination_property, $self->destination->name
                ));
            }
        }
    }
}


1;
