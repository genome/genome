package Genome::WorkflowBuilder::Event;

use strict;
use warnings;

use Genome;


class Genome::WorkflowBuilder::Event {
    is => 'Genome::WorkflowBuilder::Detail::Operation',

    has => [
        event_id => {
            is => 'String',
        },
        input_properties => {
            is => 'Text',
            is_many => 1,
        },
        output_properties => {
            is => 'Text',
            is_many => 1,
        },
    ],
    has_optional => [
        lsf_queue => {
            is => 'String',
        },
        lsf_project => {
            is => 'String',
        },
        lsf_resource => {
            is => 'String',
        },
    ],
};

# ------------------------------------------------------------------------------
# Inherited Methods
# ------------------------------------------------------------------------------

sub from_xml_element {
    my ($class, $element) = @_;

    my @input_properties = map {$_->textContent}
        $element->findnodes('.//inputproperty');
    my @output_properties = map {$_->textContent}
        $element->findnodes('.//outputproperty');

    return $class->create(
        input_properties => \@input_properties,
        output_properties => \@output_properties,
        name => $element->getAttribute('name'),
        $class->operationtype_attributes_from_xml_element($element),
    );
}

sub expected_attributes {
    return (
        event_id => 'eventId',
        lsf_project => 'lsfProject',
        lsf_queue => 'lsfQueue',
        lsf_resource => 'lsfResource',
    );
}

sub is_input_property {
    my ($self, $name) = @_;

    return Set::Scalar->new($self->input_properties)->contains($name);
}

sub is_output_property {
    my ($self, $name) = @_;

    return Set::Scalar->new($self->output_properties)->contains($name);
}

sub is_many_property {}

sub notify_input_link {
    my ($self, $link) = @_;

    unless ($self->is_input_property($link->destination_property)) {
        $self->add_input_property($link->destination_property);
    }

    return;
}

sub notify_output_link {
    my ($self, $link) = @_;

    unless ($self->is_output_property($link->source_property)) {
        $self->add_output_property($link->source_property);
    }

    return;
}

sub operation_type_attributes {
    my $self = shift;
    my %attributes;
    my %expected_attributes = $self->expected_attributes;
    for my $name (keys(%expected_attributes)) {
        my $value = $self->$name;
        if (defined($value)) {
            $attributes{$expected_attributes{$name}} = $value;
        }
    }
    return %attributes;
}

sub _execute_inline {
    my ($self, $inputs) = @_;

    die "Events are not supported for inline execution.";
}


1;
