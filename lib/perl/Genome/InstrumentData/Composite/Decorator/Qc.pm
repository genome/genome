package Genome::InstrumentData::Composite::Decorator::Qc;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Composite::Decorator::Qc {
    is => 'Genome::InstrumentData::Composite::Decorator::Base',
};

sub decorate {
    my $class = shift;
    my $operation = shift;
    my $config_name = shift;

    my $workflow = $operation->workflow_model;

    my $qc_runner_op = Workflow::Operation->create(
        name => $operation->name . ' @qc(' . $config_name . ')',
        operation_type => Workflow::OperationType::Command->get('Genome::Qc::Run'),
    );
    $qc_runner_op->workflow_model($workflow);
    my $new_input_property = 'config_name_' . $qc_runner_op->id;
    my $new_output_property = 'qc_result_' . $qc_runner_op->id;

    push @{ $workflow->operation_type->input_properties }, $new_input_property;
    push @{ $workflow->operation_type->output_properties }, $new_output_property;

    Genome::InstrumentData::Composite::Workflow::Generator::Base->_add_link_to_workflow(
        $workflow,
        left_workflow_operation_id => $workflow->get_input_connector->id,
        left_property => $new_input_property,
        right_workflow_operation_id => $qc_runner_op->id,
        right_property => 'config_name',
    );
    Genome::InstrumentData::Composite::Workflow::Generator::Base->_add_link_to_workflow(
        $workflow,
        left_workflow_operation_id => $operation->id,
        left_property => 'alignment_result',
        right_workflow_operation_id => $qc_runner_op->id,
        right_property => 'alignment_result',
    );
    Genome::InstrumentData::Composite::Workflow::Generator::Base->_add_link_to_workflow(
        $workflow,
        left_workflow_operation_id => $qc_runner_op->id,
        left_property => 'output_result',
        right_workflow_operation_id => $workflow->get_output_connector->id,
        right_property => $new_output_property,
    );

    return ('m_'. $new_input_property => $config_name);
}

1;
