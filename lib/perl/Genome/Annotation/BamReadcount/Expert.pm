package Genome::Annotation::BamReadcount::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::WorkflowBuilder::DAG;
use Genome::WorkflowBuilder::Command;
use Params::Validate qw(validate :types);

class Genome::Annotation::BamReadcount::Expert {
    is => 'Genome::Annotation::ExpertBase',
};

sub name {
    'bam-readcount';
}

sub dag {
    my $self = shift;

    my $dag = Genome::WorkflowBuilder::DAG->create(
        name => 'BamReadcount',
    );
    my $build_adaptor_op = $self->build_adaptor_op;
    $dag->add_operation($build_adaptor_op);
    $dag->connect_input(
        input_property => 'build_id',
        destination => $build_adaptor_op,
        destination_property => 'build_id',
    );
    $dag->connect_input(
        input_property => 'variant_type',
        destination => $build_adaptor_op,
        destination_property => 'variant_type',
    );

    my $run_op = $self->run_op;
    $dag->add_operation($run_op);
    $run_op->parallel_by('aligned_bam_result');
    $dag->create_link(
        source => $build_adaptor_op,
        source_property => 'bam_results',
        destination => $run_op,
        destination_property => 'aligned_bam_result',
    );
    _link(dag => $dag,
          adaptor => $build_adaptor_op,
          previous => $build_adaptor_op,
          target => $run_op,
    );

    my $annotate_op = $self->annotate_op;
    $dag->add_operation($annotate_op);
    $dag->create_link(
        source => $run_op,
        source_property => 'output_result',
        destination => $annotate_op,
        destination_property => 'readcount_results',
    );
    _link(dag => $dag,
          adaptor => $build_adaptor_op,
          previous => $build_adaptor_op,
          target => $annotate_op,
    );

    $dag->connect_output(
        output_property => 'output_result',
        source => $annotate_op,
        source_property => 'output_result',
    );

    return $dag;
}

sub build_adaptor_op {
    my $self = shift;
    return Genome::WorkflowBuilder::Command->create(
        name => 'Get inputs from build',
        command => $self->adaptor_class,
    );
}

sub run_op {
    return Genome::WorkflowBuilder::Command->create(
        name => 'Run bam-readcount',
        command => 'Genome::Annotation::BamReadcount::Run',
    );
}

sub annotate_op {
    return Genome::WorkflowBuilder::Command->create(
        name => 'Annotate vcf with readcounts',
        command => 'Genome::Annotation::BamReadcount::Annotate',
    );
}

sub _link {
    my %p = validate(@_, {
        dag => {isa => 'Genome::WorkflowBuilder::DAG'},
        adaptor => {isa => 'Genome::WorkflowBuilder::Command'},
        previous => {type => OBJECT | UNDEF},
        target => {isa => 'Genome::WorkflowBuilder::Command'},
    });

    if (defined $p{previous}) {
        $p{dag}->create_link(
            source => $p{previous},
            source_property => 'output_result',
            destination => $p{target},
            destination_property => 'input_result',
        );
    }

    for my $name ($p{target}->command->input_names) {
        next if $name eq 'input_result';
        next unless $p{adaptor}->command->can($name);
        $p{dag}->create_link(
            source => $p{adaptor},
            source_property => $name,
            destination => $p{target},
            destination_property => $name,
        );
    }
}

1;
