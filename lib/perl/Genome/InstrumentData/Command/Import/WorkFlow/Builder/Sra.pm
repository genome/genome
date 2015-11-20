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
        'retrieve source path', 'verify not imported', 'sra to bam', 'sort bam',
        'sanitize bam', 'split bam by rg', 'create instrument data',
    );
}

sub _add_sra_to_bam_op_to_workflow {
    my $self = shift;

    my $workflow = $self->_workflow;
    my $sra_to_bam_op = $self->helpers->add_operation_to_workflow_by_name($workflow, 'sra to bam');
    for my $property_mapping ( [qw/ working_directory working_directory /], [qw/ source_path sra_path /] ) {
        my ($left_property, $right_property) = @$property_mapping;
        $workflow->add_link(
            left_operation => $self->_work_flow_op_for('verify not imported'),
            left_property => $left_property,
            right_operation => $sra_to_bam_op,
            right_property => $right_property,
        );
    }

    return $sra_to_bam_op;
}

1;

