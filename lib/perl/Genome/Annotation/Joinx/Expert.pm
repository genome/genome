package Genome::Annotation::Joinx::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::WorkflowBuilder::DAG;
use Genome::WorkflowBuilder::Command;

class Genome::Annotation::Joinx::Expert {
    is => 'Genome::Annotation::ExpertBase',
};

sub name {
    'joinx';
}

sub dag {
    my $self = shift;

    my $dag = Genome::WorkflowBuilder::DAG->create(
        name => 'Joinx',
    );
    my $build_adaptor_op = $self->connected_build_adaptor_operation($dag);

    my $run_op = $self->run_op;
    $dag->add_operation($run_op);
    $self->_link(dag => $dag,
          adaptor => $build_adaptor_op,
          previous => $build_adaptor_op,
          target => $run_op,
    );

    $dag->connect_output(
        output_property => 'output_result',
        source => $run_op,
        source_property => 'output_result',
    );

    return $dag;
}

sub run_op {
    return Genome::WorkflowBuilder::Command->create(
        name => 'Run joinx',
        command => 'Genome::Annotation::Joinx::Run',
    );
}


1;
