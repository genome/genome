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

    #Make a workflow with its input and output connectors
    my $input_properties = ['instrument_data', 'reference_sequence_build', $class->_general_workflow_input_properties];
    my $tree_properties = ['name', 'params', 'version'];
    my $workflow = Workflow::Model->create(
        name => 'Speedseq',
        input_properties => [@$input_properties, @$tree_properties],
        optional_input_properties => $input_properties,
        output_properties => ['result_id'],
    );
    my $workflows = {};
    map { $workflows->{$_} = $workflow } @$alignment_objects;

    #Make a speedseq operation
    my $operation = $workflow->add_operation(
        name => 'Speedseq operation',
        operation_type => Workflow::OperationType::Command->create(
            command_class_name => 'Genome::InstrumentData::Command::AlignAndMerge',
        ),
    );

    #Connect input connectors to the operation
    my $inputs = [];
    for my $input_property (@$input_properties) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $input_property,
            right_operation => $operation,
            right_property => $input_property,
        );
        push @$inputs, ( 'm_' . $input_property => $input_data->{$input_property} );
    }
    for my $input_property (@$tree_properties) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $input_property,
            right_operation => $operation,
            right_property => $input_property,
        );
        push @$inputs, ( 'm_' . $input_property => $tree->{'action'}->[0]->{$input_property} );
    }

    #Connect output connectors to the operation
    $workflow->add_link(
        left_operation => $operation,
        left_property => 'result_id',
        right_operation => $workflow->get_output_connector,
        right_property => 'result_id',
    );

    return $workflows, $inputs;
}

# sub _workflow_

1;
