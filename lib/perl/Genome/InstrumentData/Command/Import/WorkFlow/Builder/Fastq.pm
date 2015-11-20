package Genome::InstrumentData::Command::Import::WorkFlow::Builder::Fastq;

use strict;
use warnings;

use Genome;

require List::MoreUtils;

class Genome::InstrumentData::Command::Import::WorkFlow::Builder::Fastq {
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::Builder',
};

sub _steps_to_build_workflow {
    my $self = shift;

    my @steps = ( 'retrieve source path', 'verify not imported', );
    push @steps, 'archive to fastqs' if $self->work_flow_inputs->source_files->is_tar;
    push @steps, 'fastqs to bam', 'sort bam', 'create instrument data';
 
    return @steps;
}

sub _add_archive_to_fastqs_op_to_workflow {
    my $self = shift;

    my $archive_to_fastqs_op = $self->helpers->add_operation_to_workflow_by_name($self->_workflow, 'archive to fastqs');
    return if not $archive_to_fastqs_op;
    $self->_workflow->add_link(
        left_operation => $self->_workflow->get_input_connector,
        left_property => 'working_directory',
        right_operation => $archive_to_fastqs_op,
        right_property => 'working_directory',
    );
    $self->_workflow->add_link(
        left_operation => $self->_work_flow_op_for('verify not imported'),
        left_property => 'source_path',
        right_operation => $archive_to_fastqs_op,
        right_property => 'archive_path',
    );

    return $archive_to_fastqs_op;
}

sub _add_fastqs_to_bam_op_to_workflow {
    my ($self, $previous_op) = @_;

    die 'No previous operation given!' if not $previous_op;

    my $fastqs_to_bam_op = $self->helpers->add_operation_to_workflow_by_name($self->_workflow, 'fastqs to bam');
    return if not $fastqs_to_bam_op;

    for my $property (qw/ working_directory library_name sample_name /) {
        $self->_workflow->add_link(
            left_operation => $self->_workflow->get_input_connector,
            left_property => $property,
            right_operation => $fastqs_to_bam_op,
            right_property => $property,
        );
    }
    $self->_workflow->add_link(
        left_operation => $previous_op,
        left_property => ( $previous_op->name eq 'verify not imported' ) # not ideal...
        ? 'source_path'
        : 'fastq_paths',
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
        left_property => ( $previous_op->operation_type->command_class_name->can('output_bam_path') )
        ? 'output_bam_path'
        : 'source_path',
        right_operation => $sort_bam_op,
        right_property => 'bam_path',
    );

    return $sort_bam_op;
}

1;

