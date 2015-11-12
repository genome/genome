package Genome::WorkflowBuilder::Command;

use strict;
use warnings;

use Genome;
use Cwd qw();
use Genome::Sys::LSF::ResourceParser qw(parse_lsf_params);


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

sub get_ptero_builder_task {
    require Ptero::Builder::Detail::Workflow::Task;

    my $self = shift;
    my $log_dir = shift;

    $self->validate;

    my %params = (
        name => $self->name,
        methods => [
            $self->_get_ptero_shortcut_method($log_dir),
            $self->_get_ptero_execute_method($log_dir),
        ],
    );
    if (defined $self->parallel_by) {
        $params{parallel_by} = $self->parallel_by;
    }
    return Ptero::Builder::Detail::Workflow::Task->new(%params);
}

sub _get_ptero_shortcut_method {
    require Ptero::Builder::Job;

    my $self = shift;
    my $log_dir = shift;
    return Ptero::Builder::Job->new(
        name => 'shortcut',
        service_url => Genome::Config::get('ptero_shell_command_service_url'),
        parameters => {
            commandLine => [
                'genome', 'ptero', 'wrapper',
                '--command-class', $self->command,
                '--method', 'shortcut',
                '--log-directory', $log_dir,
            ],
            environment => $self->_get_sanitized_env(),
            user => Genome::Sys->username,
            workingDirectory => Cwd::getcwd,
        },
    );
}

sub _get_ptero_execute_method {
    require Ptero::Builder::Job;


    my $self = shift;
    my $log_dir = shift;
    my $ptero_lsf_parameters = $self->_get_ptero_lsf_parameters();
    $ptero_lsf_parameters->{command} = sprintf(
        'genome ptero wrapper --command-class %s '
        .'--method execute --log-directory %s',
        $self->command, $log_dir);
    $ptero_lsf_parameters->{environment} = $self->_get_sanitized_env();
    $ptero_lsf_parameters->{user} = Genome::Sys->username;
    $ptero_lsf_parameters->{cwd} = Cwd::getcwd;

    return Ptero::Builder::Job->new(
        name => 'execute',
        service_url => Genome::Config::get('ptero_lsf_service_url'),
        parameters => $ptero_lsf_parameters);
}

sub _get_lsf_resources_from_command {
    my $self = shift;
    my $prop = $self->command->__meta__->property_meta_for_name('lsf_resource');

    if ($prop && $prop->{is_param}) {
        if ($prop->default_value) {
            return $prop->default_value;
        } else {
            die $self->command . "property lsf_resource should have a default value if it is a parameter.";
        }
    } else {
        return '';
    }
}

sub _get_ptero_lsf_parameters {
    my $self = shift;

    my $lsf_resource = $self->lsf_resource;
    if (defined($lsf_resource) && length($lsf_resource)) {
        return parse_lsf_params($lsf_resource);
    }

    return parse_lsf_params(
        $self->_get_lsf_resources_from_command
    )
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


1;
