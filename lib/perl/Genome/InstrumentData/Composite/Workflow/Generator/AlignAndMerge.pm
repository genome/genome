package Genome::InstrumentData::Composite::Workflow::Generator::AlignAndMerge;

use strict;
use warnings;
use Genome;

class Genome::InstrumentData::Composite::Workflow::Generator::AlignAndMerge {
    is => 'Genome::InstrumentData::Composite::Workflow::Generator::Base',
};

sub generate {
    my $class = shift;
    my $master_workflow = shift;
    my $block_operation = shift;
    my $tree = shift;
    my $input_data = shift;
    my $alignment_objects = shift;

    my $aligner_name = ucfirst($tree->{action}->[0]->{name});

    #Make a workflow with its input and output connectors
    my $input_properties = [$class->_general_workflow_input_properties];
    my $tree_properties = ['name', 'params', 'version'];
    my $dag_number = $class->next_counter_value;
    my $workflow = Genome::WorkflowBuilder::DAG->create(
        name => join(' ', $aligner_name, $dag_number),
    );
    my $workflows = {};
    map { $workflows->{$_} = $workflow } @$alignment_objects;

    my $command_class = 'Genome::InstrumentData::Command::AlignAndMerge';

    #Make a align_and_merge operation
    my $operation = Genome::WorkflowBuilder::Command->create(
        name => "$aligner_name operation",
        command => $command_class,
    );
    $workflow->add_operation($operation);

    my @instrument_data = map { $_->[0] } @$alignment_objects;
    my $lsf_resource_string = $command_class->lsf_resource_string_for_aligner_and_instrument_data(
        $aligner_name,
        @instrument_data
    );

    $operation->lsf_resource($lsf_resource_string);

    #Connect input connectors to the operation
    my $inputs = [];
    for my $input_property (@$input_properties) {
        $workflow->connect_input(
            input_property => $input_property,
            destination => $operation,
            destination_property => $input_property,
            is_optional => 1,
        );
        push @$inputs, ( 'm_' . $input_property => $input_data->{$input_property} );
    }
    for my $input_property (@$tree_properties) {
        my $op_input_name = $class->_operation_specific_name($workflow, $input_property);
        $workflow->connect_input(
            input_property => $op_input_name,
            destination => $operation,
            destination_property => $input_property,
            is_optional => 1,
        );
        push @$inputs, ( 'm_' . $op_input_name => $tree->{'action'}->[0]->{$input_property} );
    }

    my $reference_input_name = $tree->{'action'}->[0]->{reference};
    my $op_reference_input_name = $class->_operation_specific_name($workflow, 'reference_sequence_build');
    $workflow->connect_input(
        input_property => $op_reference_input_name,
        destination => $operation,
        destination_property => 'reference_sequence_build',
        is_optional => 1,
    );
    push @$inputs, ( 'm_' . $op_reference_input_name => $input_data->{$reference_input_name} );

    my $op_instrument_data_input_name = $class->_operation_specific_name($workflow, 'instrument_data');
    $workflow->connect_input(
        input_property => $op_instrument_data_input_name,
        destination => $operation,
        destination_property => 'instrument_data',
        is_optional => 1,
    );
    push @$inputs, ( 'm_' . $op_instrument_data_input_name => \@instrument_data );

    #Connect output connectors to the operation
    $workflow->connect_output(
        source => $operation,
        source_property => 'result_id',
        output_property => $class->_operation_specific_name($workflow, 'result_id'),
    );

    if (exists $tree->{'action'}->[0]->{decoration}) {
        push @$inputs, Genome::InstrumentData::Composite::Decorator->decorate($operation, $workflow, $tree->{'action'}->[0]->{decoration});
    }

    $class->_wire_object_workflow_to_master_workflow($master_workflow, $block_operation, $workflow);

    return $workflows, $inputs;
}

sub _wire_object_workflow_to_master_workflow {
    my $class = shift;
    my $master_workflow = shift;
    my $block_operation = shift;
    my $workflow = shift;

    $master_workflow->add_operation($workflow);

    #wire up the master to the inner workflows (just pass along the inputs and outputs)
    for my $property ($workflow->input_properties) {
        if($property eq 'force_fragment'){
            $master_workflow->create_link(
                source => $block_operation,
                source_property => $property,
                destination => $workflow,
                destination_property => $property,
            );
        }else {
            $master_workflow->connect_input(
                input_property => 'm_' . $property,
                destination => $workflow,
                destination_property => $property,
                is_optional => 1,
            );
        }
    }

    for my $property ($workflow->output_properties) {
            $master_workflow->connect_output(
            source => $workflow,
            source_property => $property,
            output_property => 'm_' . $property,
        );
    }

    return 1;
}

sub _operation_specific_name {
    my $class = shift;
    my $op = shift;
    my $name = shift;

    return join('-', $name, $op->name);
}

1;
