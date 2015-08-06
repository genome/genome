package Genome::InstrumentData::Composite::Workflow::Generator::AlignAndMerge;

use strict;
use warnings;
use Genome;

class Genome::InstrumentData::Composite::Workflow::Generator::AlignAndMerge {
    is => 'Genome::InstrumentData::Composite::Workflow::Generator::Base',
};

sub generate {
    my $class = shift;
    my $tree = shift;
    my $input_data = shift;
    my $alignment_objects = shift;

    my $aligner_name = ucfirst($tree->{action}->[0]->{name});

    #Make a workflow with its input and output connectors
    my $input_properties = ['instrument_data', 'reference_sequence_build', $class->_general_workflow_input_properties];
    my $tree_properties = ['name', 'params', 'version'];
    my $workflow = Workflow::Model->create(
        name => $aligner_name,
        input_properties => [@$input_properties, @$tree_properties],
        optional_input_properties => $input_properties,
        output_properties => ['result_id'],
    );
    my $workflows = {};
    map { $workflows->{$_} = $workflow } @$alignment_objects;

    my $command_class = 'Genome::InstrumentData::Command::AlignAndMerge';

    #Make a align_and_merge operation
    my $operation = $workflow->add_operation(
        name => "$aligner_name operation",
        operation_type => Workflow::OperationType::Command->create(
            command_class_name => $command_class,
        ),
    );

    my $instrument_data = $input_data->{instrument_data};
    my $lsf_resource_string = $command_class->lsf_resource_string_for_aligner_and_instrument_data(
        $aligner_name,
        @$instrument_data
    );

    $operation->operation_type->lsf_resource($lsf_resource_string);

    #Connect input connectors to the operation
    my $inputs = [];
    for my $input_property (@$input_properties) {
        $class->_add_link_to_workflow(
            $workflow,
            left_workflow_operation_id => $workflow->get_input_connector->id,
            left_property => $input_property,
            right_workflow_operation_id => $operation->id,
            right_property => $input_property,
        );
        push @$inputs, ( 'm_' . $input_property => $input_data->{$input_property} );
    }
    for my $input_property (@$tree_properties) {
        $class->_add_link_to_workflow(
            $workflow,
            left_workflow_operation_id => $workflow->get_input_connector->id,
            left_property => $input_property,
            right_workflow_operation_id => $operation->id,
            right_property => $input_property,
        );
        push @$inputs, ( 'm_' . $input_property => $tree->{'action'}->[0]->{$input_property} );
    }

    #Connect output connectors to the operation
    $class->_add_link_to_workflow(
        $workflow,
        left_workflow_operation_id => $operation->id,
        left_property => 'result_id',
        right_workflow_operation_id => $workflow->get_output_connector->id,
        right_property => 'result_id',
    );

    if (exists $tree->{'action'}->[0]->{decoration}) {
        push @$inputs, Genome::InstrumentData::Composite::Decorator->decorate($operation, $tree->{'action'}->[0]->{decoration});
    }

    return $workflows, $inputs;
}

1;
