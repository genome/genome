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

    my $name = $operation->name . ' @qc(' . $config_name . ')';
    my $qc_runner_op = $class->create_qc_runner_op($name);
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
        left_workflow_operation_id => $workflow->get_input_connector->id,
        left_property => 'result_users',
        right_workflow_operation_id => $qc_runner_op->id,
        right_property => 'result_users',
    );
    $class->add_alignment_link($workflow, $operation, $qc_runner_op);
    Genome::InstrumentData::Composite::Workflow::Generator::Base->_add_link_to_workflow(
        $workflow,
        left_workflow_operation_id => $qc_runner_op->id,
        left_property => 'output_result',
        right_workflow_operation_id => $workflow->get_output_connector->id,
        right_property => $new_output_property,
    );

    return ('m_'. $new_input_property => $config_name);
}

sub create_qc_runner_op {
    my $class = shift;
    my $name = shift;

    my $qc_runner_op = Workflow::Operation->create(
        name => $name,
        operation_type => Workflow::OperationType::Command->get('Genome::Qc::Run'),
    );

    return $qc_runner_op;
}

sub add_alignment_link {
    my ($class, $workflow, $operation, $qc_runner_op) = @_;

    Genome::InstrumentData::Composite::Workflow::Generator::Base->_add_link_to_workflow(
        $workflow,
        left_workflow_operation_id => $operation->id,
        left_property => 'alignment_result',
        right_workflow_operation_id => $qc_runner_op->id,
        right_property => 'alignment_result',
    );
}

1;
