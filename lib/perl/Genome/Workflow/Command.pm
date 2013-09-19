package Genome::Workflow::Command;

use strict;
use warnings;

use Genome;


class Genome::Workflow::Command {
    is => 'Genome::Workflow::Detail::Operation',

    has => [
        command => {
            is => 'Command',
        },
    ],
};


# ------------------------------------------------------------------------------
# Inherited Methods
# ------------------------------------------------------------------------------

sub input_properties {
    my $self = shift;
    return map {$_->property_name} $self->command->__meta__->properties(
        is_input => 1, is_optional => 0);
}

my %_EXPECTED_ATTRIBUTES = (
    lsf_project => 'lsfProject',
    lsf_queue => 'lsfQueue',
    lsf_resource => 'lsfResource',
);
sub operation_type_attributes {
    my $self = shift;
    my %attributes = (
        commandClass => $self->command,
    );
    for my $command_property (keys(%_EXPECTED_ATTRIBUTES)) {
        my $value = $self->_get_attribue_from_command($command_property);
        if (defined($value)) {
            $attributes{$_EXPECTED_ATTRIBUTES{$command_property}} = $value;
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


1;
