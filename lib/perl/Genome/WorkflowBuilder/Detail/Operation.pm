package Genome::WorkflowBuilder::Detail::Operation;

use strict;
use warnings;

use Genome;
use Genome::WorkflowBuilder::Detail::TypeMap;
use IO::Scalar qw();
use Set::Scalar qw();
use XML::LibXML qw();
use Carp qw(confess);
use Data::Dumper qw();
use List::MoreUtils qw(firstval);
use Params::Validate qw(validate_pos :types);
use Path::Class qw();

class Genome::WorkflowBuilder::Detail::Operation {
    is => 'Genome::WorkflowBuilder::Detail::Element',
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
    has_transient => {
        constant_values => {
            is => 'HASH',
            default => {},
        },
        outputs => {
            is => 'HASH',
            default => {},
        }
    }
};


# ------------------------------------------------------------------------------
# Abstract methods
# ------------------------------------------------------------------------------

sub from_xml_element {
    my ($class, $element, @rest) = @_;

    # Prevent accidental recursion when subclasses don't override this method
    unless ($class eq 'Genome::WorkflowBuilder::Detail::Operation') {
        confess $class->error_message(sprintf(
                "from_xml_element not implemented in subclass %s", $class));
    }

    my $has_file = $element->hasAttribute('workflowFile');
    if ($has_file) {
        my $file = $class->_find_inner_workflow_file_for_element($element);
        my $result = $class->from_xml_filename($file);
        if ($element->hasAttribute('name')) {
            $result->name($element->getAttribute('name'));
        }
        return $result;
    } else {
        my $subclass = $class->_get_subclass_from_element($element);
        return $subclass->from_xml_element($element, @rest);
    }
}

sub input_properties {}
sub optional_input_properties {}
sub output_properties {}
sub operation_type_attributes {}

sub is_input_property {}
sub is_optional_input_property {}
sub is_output_property {}
sub is_many_property {}


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
    my ($class, $filename) = @_;

    my $doc = XML::LibXML->load_xml(location => $filename);
    return $class->from_xml_element($doc->documentElement);
}

sub operation_type {
    my $self = shift;

    return Genome::WorkflowBuilder::Detail::TypeMap::type_from_class($self->class);
}

sub declare_constant {
    my $self = shift;
    my %constants = @_;

    while (my ($key, $value) = each %constants) {
        unless ($self->is_input_property($key)) {
            die sprintf("No input named (%s) on operation named (%s)",
                $key, $self->name);
        }
        $self->constant_values->{$key} = $value;
    }
}

# ------------------------------------------------------------------------------
# Inherited methods
# ------------------------------------------------------------------------------

sub notify_input_link {}

sub notify_output_link {}

sub get_xml {
    my $self = shift;

    $self->validate;

    my $doc = XML::LibXML::Document->new();
    $doc->setDocumentElement($self->get_xml_element);

    return $doc->toString(1);
}

sub remove_links_from_unspecified_input {
    return;
}

sub remove_links_to_unused_input_properties {
    return;
}

sub recurse_do {
    return;
}

sub get_xml_element {
    my $self = shift;

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
        die sprintf("Operation name '%s' is not allowed", $self->name);
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

    $self->_add_property_xml_elements($element);

    my %attributes = $self->operation_type_attributes;
    for my $attr_name (keys(%attributes)) {
        $element->setAttribute($attr_name, $attributes{$attr_name});
    }

    return $element;
}

sub _add_property_xml_elements {
    my ($self, $element) = @_;

    my @input_properties = sort Set::Scalar->new($self->input_properties,
        $self->optional_input_properties)->members;
    map {$self->_add_property_xml_element($element, 'inputproperty', $_)}
        @input_properties;
    map {$self->_add_property_xml_element($element, 'outputproperty', $_)}
        $self->output_properties;

    return;
}

sub _add_property_xml_element {
    my ($self, $element, $xml_tag, $text) = @_;

    my $inner_element = XML::LibXML::Element->new($xml_tag);
    $inner_element->appendText($text);
    if ($self->is_optional_input_property($text)) {
        $inner_element->setAttribute('isOptional', 'Y');
    }
    $element->addChild($inner_element);

    return;
}

sub _find_inner_workflow_file_for_element {
    my ($class, $element) = @_;

    my $file = $element->getAttribute('workflowFile');
    unless ($file) {
        $class->fatal_message('No workflow file set on element.');
    }

    my $path = Path::Class::file($file);

    unless ($path->is_absolute()) {
        #relative paths should be relative to the current WF file, rather than the CWD.
        my $current_file = $element->ownerDocument->URI;
        if (defined $current_file and -e $current_file) {
            my $current_path = Path::Class::file($current_file);
            $path = $current_path->dir->file($file);
            $file = $path->stringify();
        } else {
            $class->fatal_message('Cannot resolve relative path <%s> in inner workflow.', $file);
        }
    }

    unless (-e $file) {
        $class->fatal_message('Inner workflow file <%s> not found!', $file);
    }

    return $file;
}

sub _get_subclass_from_element {
    my ($class, $element) = @_;
    my $nodes = $element->find('operationtype');
    my $operation_type_element = $nodes->pop;
    return Genome::WorkflowBuilder::Detail::TypeMap::class_from_type(
        $operation_type_element->getAttribute('typeClass'));
}

sub _get_sanitized_env {
    my $self = shift;

    my $env = {};
    while (my($key, $value) = each %ENV) {
        if (defined($value)) {
            $env->{$key} = $value;
        } else {
            $self->warning_message("Found Environment variable with no value: %s", $key);
        }
    }
    return $env;
}

sub operationtype_attributes_from_xml_element {
    my ($class, $element) = @_;

    my %properties;
    my %expected_attributes = $class->expected_attributes;
    for my $property_name (keys(%expected_attributes)) {
        my $attribute_name = $expected_attributes{$property_name};
        $properties{$property_name} = $class->_get_value_from_xml_element(
            $element, $attribute_name);
    }
    return %properties;
}

sub _get_value_from_xml_element {
    my ($class, $element, $name) = @_;

    my $nodes = $element->find('operationtype');
    my $operation_type_element = $nodes->pop;
    return $operation_type_element->getAttribute($name);
}

sub execute_inline {
    my ($self, $inputs) = validate_pos(@_, {type => OBJECT}, {type => HASHREF});

    $self->outputs({});

    if (defined($self->parallel_by)) {
        for my $parallel_input (@{$inputs->{$self->parallel_by}}) {
            my %loop_inputs = %{$inputs};
            $loop_inputs{$self->parallel_by} = $parallel_input;
            my $loop_outputs = $self->_execute_inline(\%loop_inputs);
            $self->_store_outputs($loop_outputs);
        }
    } else {
        my $outputs = $self->_execute_inline($inputs);
        $self->_store_outputs($outputs);
    }

    return $self->outputs;
}

sub _execute_inline {
    my $self = shift;
    my $inputs = shift;

    die "Abstract";
}

sub _store_outputs {
    my ($self, $outputs) = validate_pos(@_, {type => OBJECT}, {type => HASHREF});

    if (defined($self->parallel_by)) {
        my $self_outputs = $self->outputs;
        while (my ($key, $value) = each %{$outputs}) {
            push @{$self_outputs->{$key}}, $value;
        }
    } else {
        $self->outputs($outputs);
    }
    return;
}

1;
