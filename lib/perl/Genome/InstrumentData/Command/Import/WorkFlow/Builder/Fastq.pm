package Genome::InstrumentData::Command::Import::WorkFlow::Builder::Fastq;

use strict;
use warnings;

use Genome;

require List::MoreUtils;

class Genome::InstrumentData::Command::Import::WorkFlow::Builder::Fastq {
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::Builder',
};

sub _steps_to_build_workflow {
    return ( 'verify not imported', 'fastqs to bam', 'create instrument data' );
}

sub _add_fastqs_to_bam_op_to_workflow {
    my ($self, $previous_op) = @_;

    die 'No previous operation given!' if not $previous_op;

    my $name = 'fastqs to bam';
    my $fastqs_to_bam_op = Genome::WorkflowBuilder::Command->create(
        name => $name,
        command => $self->work_flow_operation_class_for_name($name),
    );
    $self->_dag->add_operation($fastqs_to_bam_op);

    for my $property (qw/ working_directory library /) {
        $self->_dag->connect_input(
            input_property => $property,
            destination => $fastqs_to_bam_op,
            destination_property => $property,
        );
    }

    $self->_dag->create_link(
        source => $previous_op,
        source_property => 'output_paths',
        destination => $fastqs_to_bam_op,
        destination_property => 'fastq_paths',
    );

    return $fastqs_to_bam_op;
}

1;

