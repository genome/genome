package Genome::WorkflowBuilder::Command;

use strict;
use warnings;

use Genome;
use Cwd qw();
use Genome::Sys::LSF::ResourceParser qw(parse_lsf_params);
use Carp qw();
use Data::Dump qw(pp);
use Try::Tiny;
use File::Spec;


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
        Carp::confess(sprintf("Failed to load command class (%s): %s\n",
                $self->command, $error));
    }
    return $self;
}

# ------------------------------------------------------------------------------
# Inherited Methods
# ------------------------------------------------------------------------------

sub from_xml_element {
    my ($class, $element) = @_;

    return $class->create(
        name => $element->getAttribute('name'),
        parallel_by => $element->getAttribute('parallelBy'),
        $class->operationtype_attributes_from_xml_element($element),
    );
}

sub expected_attributes {
    return (
        command => 'commandClass',
        lsf_project => 'lsfProject',
        lsf_queue => 'lsfQueue',
        lsf_resource => 'lsfResource',
    );
}

sub input_properties {
    my $self = shift;

    my @metas = $self->command->__meta__->properties(
        is_input => 1, is_optional => 0);

    my @metas_without_defaults = grep {! defined($_->default_value)} @metas;

    my @result = map {$_->property_name} @metas_without_defaults;
    return sort @result;
}

sub operation_type_attributes {
    my $self = shift;
    my %attributes = (
        commandClass => $self->command,
    );
    my %expected_attributes = $self->expected_attributes;
    for my $name (keys(%expected_attributes)) {
        my $value;
        if (defined($self->$name)) {
            $value = $self->$name;
        } else {
            $value = $self->_get_attribute_from_command($name);
        }

        if (defined($value)) {
            $attributes{$expected_attributes{$name}} = $value;
        }
    }

    unless (defined $attributes{'lsfQueue'}) {
        $attributes{'lsfQueue'} = Genome::Config::get('lsf_queue_build_worker');
    }

    return %attributes;
}

sub output_properties {
    my $self = shift;
    return sort map {$_->property_name} $self->command->__meta__->properties(
        is_output => 1);
}

sub validate {
    my $self = shift;

    $self->SUPER::validate();

    if (defined($self->parallel_by)) {
        my $prop = $self->command->__meta__->properties(
            property_name => $self->parallel_by, is_input => 1);
        if (!defined($prop)) {
            die sprintf("Failed to verify that requested " .
                    "parallel_by property '%s' was an input " .
                    "on command (%s)",
                    $self->parallel_by, $self->command->class);
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

sub _get_attribute_from_command {
    my ($self, $property_name) = @_;

    my $property = $self->command->__meta__->properties(
        property_name => $property_name);
    return unless defined $property;

    if (defined $property->default_value) {
        return $property->default_value;
    }
    elsif ($property->calculated_default) {
        return $property->calculated_default->();
    }
    else {
        return;
    }
}

sub _execute_inline {
    my ($self, $inputs) = @_;

    my $cmd = $self->_instantiate_command($inputs);
    $self->_run_command($cmd);
    return _get_command_outputs($cmd, $self->command);
}

sub _instantiate_command {
    my ($self, $inputs) = @_;

    $self->status_message("Instantiating command %s", $self->command);

    my $pkg = $self->command;
    my $cmd = try {
        eval "use $pkg";
        $pkg->create(%$inputs)
    } catch {
        Carp::confess sprintf(
            "Failed to instantiate class (%s) with inputs (%s): %s",
            $pkg, pp($inputs), pp($_))
    };

    return $cmd;
}

sub _run_command {
    my ($self, $cmd) = @_;

    $self->status_message("Running command %s", $self->command);

    my $ret = try {
        $cmd->execute()
    } catch {
        Carp::confess sprintf(
            "Crashed in execute for command %s: %s",
            $self->command, $_,
        );
    };
    unless ($ret) {
        Carp::confess sprintf("Failed to execute for command %s.",
            $self->command,
        );
    }

    $self->status_message("Succeeded to execute command %s", $self->command);
}

sub _get_command_outputs {
    my ($cmd, $pkg) = @_;

    my %outputs;
    my @output_properties = $pkg->__meta__->properties(is_output => 1);
    for my $prop (@output_properties) {
        my $prop_name = $prop->property_name;
        my $value = $prop->is_many ? [$cmd->$prop_name] : $cmd->$prop_name;
        $outputs{$prop_name} = $value;
    }

    return \%outputs;
}


1;
