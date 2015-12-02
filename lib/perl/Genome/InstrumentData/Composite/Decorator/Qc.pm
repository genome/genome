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
    my $workflow = shift;
    my $config_name = shift;

    my $name = $operation->name . ' @qc(' . $config_name . ')';
    my $qc_runner_op = $class->create_qc_runner_op($name);
    $workflow->add_operation($qc_runner_op);
    my $new_input_property = 'config_name_' . $qc_runner_op->id;
    my $new_output_property = 'qc_result_' . $qc_runner_op->id;

    $workflow->connect_input(
        input_property => $new_input_property,
        destination => $qc_runner_op,
        destination_property => 'config_name',
    );
    $workflow->connect_input(
        input_property => 'result_users',
        destination => $qc_runner_op,
        destination_property => 'result_users',
    );
    $class->add_alignment_link($workflow, $operation, $qc_runner_op);
    $workflow->connect_output(
        source => $qc_runner_op,
        source_property => 'output_result',
        output_property => $new_output_property,
    );

    return ('m_'. $new_input_property => $config_name);
}

sub create_qc_runner_op {
    my $class = shift;
    my $name = shift;

    my $qc_runner_op = Genome::WorkflowBuilder::Command->create(
        name => $name,
        command => 'Genome::Qc::Run',
    );

    return $qc_runner_op;
}

sub add_alignment_link {
    my ($class, $workflow, $operation, $qc_runner_op) = @_;

    $workflow->create_link(
        source => $operation,
        source_property => 'alignment_result',
        destination => $qc_runner_op,
        destination_property => 'alignment_result',
    );
}

1;
