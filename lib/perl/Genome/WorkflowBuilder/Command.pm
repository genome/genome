package Genome::WorkflowBuilder::Command;

use strict;
use warnings;

use Genome;


class Genome::WorkflowBuilder::Command {
    is => 'Genome::WorkflowBuilder::Detail::Operation',

    has => [
        command => {
            is => 'Command',
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

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    eval sprintf("require %s", $self->command);
    my $error = $@;
    if ($error) {
        Carp::confess(sprintf("Failed to load command class (%s)",
                $self->command));
    }
    return $self;
}


# ------------------------------------------------------------------------------
# Inherited Methods
# ------------------------------------------------------------------------------

sub from_xml_element {
    my ($class, $element) = @_;

    my $command_class = $class->_get_command_class_from_xml_element($element);
    return $class->create(
        name => $element->getAttribute('name'),
        command => $command_class,
    );
}

my %_EXPECTED_ATTRIBUTES = (
    lsf_project => 'lsfProject',
    lsf_queue => 'lsfQueue',
    lsf_resource => 'lsfResource',
);
sub input_properties {
    my $self = shift;
    my @result = map {$_->property_name} $self->command->__meta__->properties(
        is_input => 1, is_optional => 0);
    push @result, grep {!exists $_EXPECTED_ATTRIBUTES{$_}} map {$_->property_name} $self->command->__meta__->properties(
        is_param => 1, is_optional => 0);
    return @result;
}

sub operation_type_attributes {
    my $self = shift;
    my %attributes = (
        commandClass => $self->command,
    );
    for my $name (keys(%_EXPECTED_ATTRIBUTES)) {
        my $value;
        if (defined($self->$name)) {
            $value = $self->$name;
        } else {
            $value = $self->_get_attribue_from_command($name);
        }

        if (defined($value)) {
            $attributes{$_EXPECTED_ATTRIBUTES{$name}} = $value;
        }
    }
    return %attributes;
}

sub output_properties {
    my $self = shift;
    return map {$_->property_name} $self->command->__meta__->properties(
        is_output => 1);
}

sub validate {
    my $self = shift;

    $self->SUPER::validate();

    if (defined($self->parallel_by)) {
        my $prop = $self->command->__meta__->properties(
            property_name => $self->parallel_by, is_input => 1);
        if (!defined($prop)) {
            die $self->error_message(sprintf("Failed to verify that requested "
                    . "parallel_by property '%s' was an input",
                    $self->parallel_by));
        }
    }
}

sub is_input_property {
    my ($self, $property_name) = @_;
    return $self->command->__meta__->properties(
            property_name => $property_name, is_input => 1)
        || $self->command->__meta__->properties(
            property_name => $property_name, is_param => 1);
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


# ------------------------------------------------------------------------------
# Private Methods
# ------------------------------------------------------------------------------

sub _get_attribue_from_command {
    my ($self, $property_name) = @_;

    my $property = $self->command->__meta__->properties(
        property_name => $property_name);
    if (defined($property)) {
        return $property->default_value;
    } else {
        return;
    }
}

sub _get_command_class_from_xml_element {
    my ($class, $element) = @_;

    my $nodes = $element->find('operationtype');
    my $operation_type_element = $nodes->pop;
    return $operation_type_element->getAttribute('commandClass');
}


1;
