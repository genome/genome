package Genome::InstrumentData::Composite::Decorator::AlignAndMergeQc;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Composite::Decorator::AlignAndMergeQc {
    is => 'Genome::InstrumentData::Composite::Decorator::Qc',
};

sub create_qc_runner_op {
    my $class = shift;

    my $qc_runner_op = $class->SUPER::create_qc_runner_op(@_);
    $qc_runner_op->parallel_by('alignment_result');

    return $qc_runner_op;
}

sub add_alignment_link {
    my $class = shift;
    my $workflow = shift;
    my $operation = shift;
    my $qc_runner_op = shift;

    $workflow->create_link(
        source => $operation,
        source_property => 'per_lane_alignment_results',
        destination => $qc_runner_op,
        destination_property => 'alignment_result',
    );
}

1;
