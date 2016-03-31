package Genome::InstrumentData::Command::Import::WorkFlow::Builder::Sra;

use strict;
use warnings;

use Genome;

require List::MoreUtils;

class Genome::InstrumentData::Command::Import::WorkFlow::Builder::Sra {
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::Builder',
};

sub _steps_to_build_workflow {
    return (
        'verify not imported', 'sra to bam', 'sort bam',
        'sanitize and split bam', 'create instrument data',
    );
}

sub _add_sra_to_bam_op_to_workflow {
    my ($self, $previous_op) = @_;

    die 'No previous op given to _add_sra_to_bam_op_to_workflow!' if not $previous_op;

    my $name = 'sra to bam';
    my $sra_to_bam_op = Genome::WorkflowBuilder::Command->create(
        name => $name,
        command => $self->work_flow_operation_class_for_name($name),
    );
    $self->_dag->add_operation($sra_to_bam_op);

    $self->_dag->connect_input(
        input_property => 'working_directory',
        destination => $sra_to_bam_op,
        destination_property => 'working_directory',
    );
    $self->_dag->create_link(
        source => $previous_op,
        source_property => 'output_path',
        destination => $sra_to_bam_op,
        destination_property => 'sra_path',
    );

    return $sra_to_bam_op;
}

1;

