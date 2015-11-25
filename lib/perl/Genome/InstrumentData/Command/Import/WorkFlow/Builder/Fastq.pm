package Genome::InstrumentData::Command::Import::WorkFlow::Builder::Fastq;

use strict;
use warnings;

use Genome;

require List::MoreUtils;

class Genome::InstrumentData::Command::Import::WorkFlow::Builder::Fastq {
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::Builder',
};

sub _steps_to_build_workflow {
    return ( 'retrieve source path', 'verify not imported', 'fastqs to bam', 'sort bam', 'create instrument data' );
}

sub _add_fastqs_to_bam_op_to_workflow {
    my ($self, $previous_op) = @_;

    die 'No previous operation given!' if not $previous_op;

    my $fastqs_to_bam_op = $self->helpers->add_operation_to_workflow_by_name($self->_workflow, 'fastqs to bam');
    return if not $fastqs_to_bam_op;

    for my $property (qw/ working_directory library /) {
        $self->_workflow->add_link(
            left_operation => $self->_workflow->get_input_connector,
            left_property => $property,
            right_operation => $fastqs_to_bam_op,
            right_property => $property,
        );
    }
    $self->_workflow->add_link(
        left_operation => $previous_op,
        left_property => 'source_path',
        right_operation => $fastqs_to_bam_op,
        right_property => 'fastq_paths',
    );

    return $fastqs_to_bam_op;
}

sub _add_sort_bam_op_to_workflow {
    my ($self, $previous_op) = @_;

    die 'No previous op given to _add_sort_bam_op_to_workflow!' if not $previous_op;

    my $sort_bam_op = $self->helpers->add_operation_to_workflow_by_name($self->_workflow, 'sort bam');
    $self->_workflow->add_link(
        left_operation => $previous_op,
        left_property => 'output_bam_path',
        right_operation => $sort_bam_op,
        right_property => 'bam_path',
    );

    return $sort_bam_op;
}

1;

