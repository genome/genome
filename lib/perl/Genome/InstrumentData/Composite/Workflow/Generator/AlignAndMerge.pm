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
    my $input_properties = ['instrument_data', 'reference_sequence_build', $class->_general_workflow_input_properties];
    my $tree_properties = ['name', 'params', 'version'];
    my $workflow = Genome::WorkflowBuilder::DAG->create(
        name => $aligner_name,
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

    my $instrument_data = $input_data->{instrument_data};
    my $lsf_resource_string = $command_class->lsf_resource_string_for_aligner_and_instrument_data(
        $aligner_name,
        @$instrument_data
    );

    $operation->lsf_resource($lsf_resource_string);

    #Connect input connectors to the operation
    my $inputs = [];
    for my $input_property (@$input_properties) {
        $class->_add_link_to_workflow(
            $workflow,
            source => $workflow,
            source_property => $input_property,
            destination => $operation,
            destination_property => $input_property,
        );
        push @$inputs, ( 'm_' . $input_property => $input_data->{$input_property} );
    }
    for my $input_property (@$tree_properties) {
        $class->_add_link_to_workflow(
            $workflow,
            source => $workflow,
            source_property => $input_property,
            destination => $operation,
            destination_property => $input_property,
        );
        push @$inputs, ( 'm_' . $input_property => $tree->{'action'}->[0]->{$input_property} );
    }

    #Connect output connectors to the operation
    $class->_add_link_to_workflow(
        $workflow,
        source => $operation,
        source_property => 'result_id',
        destination => $workflow,
        destination_property => 'result_id',
    );

    if (exists $tree->{'action'}->[0]->{decoration}) {
        push @$inputs, Genome::InstrumentData::Composite::Decorator->decorate($operation, $tree->{'action'}->[0]->{decoration});
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
            $class->_add_link_to_workflow($master_workflow,
                source => $block_operation,
                source_property => $property,
                destination => $workflow,
                destination_property => $property,
            );
        }else {
            $class->_add_link_to_workflow($master_workflow,
                source => $master_workflow,
                source_property => 'm_' . $property,
                destination => $workflow,
                destination_property => $property,
            );
        }
    }

    for my $property ($workflow->output_properties) {
        $class->_add_link_to_workflow($master_workflow,
            source => $workflow,
            source_property => $property,
            destination => $master_workflow,
            destination_property => 'm_' . $property,
        );
    }

    return 1;
}

1;
