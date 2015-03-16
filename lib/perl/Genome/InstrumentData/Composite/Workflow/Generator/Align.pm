package Genome::InstrumentData::Composite::Workflow::Generator::Align;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Composite::Workflow::Generator::Align {
    is => 'Genome::InstrumentData::Composite::Workflow::Generator::Base',
};

sub generate {
    my $class = shift;
    my $tree = shift;
    my $input_data = shift;
    my $alignment_objects = shift;

    my $workflows = {};
    my $inputs = [];

    for my $obj (@$alignment_objects) {
        my ($workflow, $input) = $class->_generate_workflow_for_instrument_data(
            $tree, $input_data, @$obj,
        );
        $workflows->{$obj} = $workflow;
        push @$inputs, @$input;
    }

    return $workflows, $inputs;
}

sub _generate_workflow_for_instrument_data {
    my $class = shift;
    my $tree = shift;
    my $input_data = shift;
    my $instrument_data = shift;
    my %options = @_; #segment/read selection

    #First, get all the individual operations and figure out all the inputs/outputs we'll need
    my @operations = ();
    for my $subtree (@{ $tree->{action} } ) {
        push @operations, $class->_create_operations_for_alignment_tree($subtree, $instrument_data, %options);
    }

    my @input_properties = (
        $class->_general_workflow_input_properties(),
        $class->_instrument_data_workflow_input_properties($instrument_data, %options),
        (map { $class->_input_properties_for_operation($_) } @operations)
    );

    my @output_properties;
    for my $leaf (@{ $tree->{action} }) {
        push @output_properties, join('_', 'result_id', $leaf->{$class->_operation_key($instrument_data, %options)}->name);
    }

    #Next create the model, and add all the operations to it
    my $workflow_name = 'Alignment Dispatcher for ' . $instrument_data->id;
    if(exists $options{instrument_data_segment_id}) {
        $workflow_name .= ' (segment ' . $options{instrument_data_segment_id} . ')';
    }

    my $workflow = Workflow::Model->create(
        name => $workflow_name,
        input_properties => \@input_properties,
        optional_input_properties => \@input_properties,
        output_properties => \@output_properties,
    );

    for my $op (@operations) {
        $op->workflow_model($workflow);
    }

    #Last wire up the workflow
    my $inputs_for_links = $class->_generate_alignment_workflow_links($workflow, $tree, $input_data, $instrument_data, %options);

    return $workflow, $inputs_for_links;
}

sub _create_operations_for_alignment_tree {
    my $class = shift;
    my $tree = shift;
    my $instrument_data = shift;
    my %options = @_;

    my @operations = ();

    return @operations if exists $tree->{$class->_operation_key($instrument_data, %options)}; #generation already complete for this subtree

    push @operations, $class->_generate_operation($tree, $instrument_data, %options);
    if(exists $tree->{parent}) {
        push @operations, $class->_create_operations_for_alignment_tree($tree->{parent}, $instrument_data, %options);
    }

    return @operations;
}

sub _input_properties_for_operation {
    my $class = shift;
    my $operation = shift;

    my @input_properties = ();
    my $op_name = $operation->name;
    my $op_type = $operation->operation_type;

    for my $prop ('reference_build_id', 'annotation_build_id') {
        if(grep $_ eq $prop, @{ $op_type->input_properties }) {
            push @input_properties, join('_', $prop, $op_name);
        }
    }

    push @input_properties, join('_', 'name', $op_name), join('_', 'params', $op_name), join('_', 'version', $op_name);

    return @input_properties;
}

sub _generate_alignment_workflow_links {
    my $class = shift;
    my $workflow = shift;
    my $tree = shift;
    my $input_data = shift;
    my $instrument_data = shift;
    my %options = @_;

    my $inputs = [];
    for my $subtree (@{ $tree->{action} }) {
        my $input = $class->_create_links_for_subtree($workflow, $subtree, $input_data, $instrument_data, %options);
        push @$inputs, @$input;

        #create link for final result of that subtree
        my $op = $subtree->{$class->_operation_key($instrument_data, %options)};
        $class->_add_link_to_workflow($workflow,
            left_workflow_operation_id => $op->id,
            left_property => 'result_id',
            right_workflow_operation_id => $workflow->get_output_connector->id,
            right_property => join('_', 'result_id', $op->name),
        );
    }

    return $inputs;
}

sub _create_links_for_subtree {
    my $class = shift;
    my $workflow = shift;
    my $tree = shift;
    my $input_data = shift;
    my $instrument_data = shift;
    my %options = @_;

    my $link_key = $class->_operation_key($instrument_data, %options) . '_links';
    return 1 if exists $tree->{$link_key};
    $tree->{$link_key} = 1;

    my $operation = $tree->{$class->_operation_key($instrument_data, %options)};
    die('No operation on tree!') unless $operation;

    my @properties = ('name', 'params', 'version');
    push @properties, 'reference_build_id'
        if ($tree->{type} eq 'align');
    push @properties, 'annotation_build_id'
        if (defined $tree->{annotation});

    my $inputs = [];
    for my $property (@properties) {
        my $left_property = join('_', $property, $operation->name);
        $class->_add_link_to_workflow($workflow,
            left_workflow_operation_id => $workflow->get_input_connector->id,
            left_property => $left_property,
            right_workflow_operation_id => $operation->id,
            right_property => $property,
        );

        my $value;

        if($property eq 'reference_build_id') {
            my $reference = $input_data->{$tree->{reference}};
            $value = $reference->id;
            unless($value) {
                die $class->error_message('Found no reference for ' . $tree->{reference});
            }
        } elsif ($property eq 'annotation_build_id') {
            my $annotation_build = $input_data->{$tree->{annotation}};
            $value = $annotation_build->id;
            unless($value) {
                die $class->error_message('Found no annotation build for ' . $tree->{annotation});
            }
        } else {
            $value = $tree->{$property};
        }

        push @$inputs, (
            'm_'.$left_property => $value,
        );
    }

    for my $property ($class->_general_workflow_input_properties) {
        $class->_add_link_to_workflow($workflow,
            left_workflow_operation_id => $workflow->get_input_connector->id,
            left_property => $property,
            right_workflow_operation_id => $operation->id,
            right_property => $property,
        );
    }

    if(exists $tree->{parent}) {
        my $parent_operation = $tree->{parent}{$class->_operation_key($instrument_data, %options)};
        $class->_add_link_to_workflow($workflow,
            left_workflow_operation_id => $parent_operation->id,
            left_property => 'output_result_id',
            right_workflow_operation_id => $operation->id,
            right_property => 'alignment_result_input_id',
        );

        return $class->_create_links_for_subtree($workflow, $tree->{parent}, $instrument_data, %options);
    } else {
        #we're at the root--so create links for passing the segment info, etc.
        for my $property ($class->_instrument_data_workflow_input_properties($instrument_data, %options)) {
            my ($simple_property_name) = $property =~ m/^(.+?)__/;
            $class->_add_link_to_workflow($workflow,
                left_workflow_operation_id => $workflow->get_input_connector->id,
                left_property => $property,
                right_workflow_operation_id => $operation->id,
                right_property => $simple_property_name,
            );

            push @$inputs, (
                'm_'.$property => $options{$simple_property_name},
            );
        }
    }

    return $inputs;
}

sub _instrument_data_workflow_input_properties {
    my $class = shift;
    my $instrument_data = shift;
    my %options = @_;

    my @input_properties = map(
        join('__', $_, $instrument_data->id, (exists $options{instrument_data_segment_id} ? $options{instrument_data_segment_id} : ())),
        ('instrument_data_id', 'instrument_data_segment_id', 'instrument_data_segment_type', 'instrument_data_filter')
    );

    return @input_properties;
}
1;
