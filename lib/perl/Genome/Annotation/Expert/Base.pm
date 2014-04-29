package Genome::Annotation::Expert::Base;

use strict;
use warnings FATAL => 'all';
use Genome;
use Params::Validate qw(validate_pos validate :types);

class Genome::Annotation::Expert::Base {
    is => 'Genome::Annotation::ComponentBase',
    is_abstract => 1,
};

sub name {
    die "Abstract";
}

sub priority {
    # Higher priority experts are run earlier.
    # Ties are broken by expert name (alphabetically early names are
    # run earlier).
    return 0;
}

sub run_class {
    my $self = shift;
    my $factory = Genome::Annotation::Factory->create();
    return $factory->get_class('runners', $self->name);
}

sub adaptor_class {
    my $self = shift;
    my $factory = Genome::Annotation::Factory->create();
    return $factory->get_class('adaptors', $self->name);
}

sub run_operation {
    my $self = shift;
    return Genome::WorkflowBuilder::Command->create(
        name => sprintf('Run %s', $self->name),
        command => $self->run_class,
    );
}

sub connected_run_operation {
    my ($self, $dag) = validate_pos(@_, 1, 1);

    my $run_operation = $self->run_operation;
    $dag->add_operation($run_operation);
    $dag->connect_input(
        input_property => 'input_result',
        destination => $run_operation,
        destination_property => 'input_result',
    );
    $dag->connect_output(
        output_property => 'output_result',
        source => $run_operation,
        source_property => 'output_result',
    );
    return $run_operation;
}

sub build_adaptor_operation {
    my $self = shift;
    return Genome::WorkflowBuilder::Command->create(
        name => 'Get inputs from build',
        command => $self->adaptor_class,
    );
}

sub connected_build_adaptor_operation {
    my ($self, $dag) = validate_pos(@_, 1, 1);

    my $build_adaptor_operation = $self->build_adaptor_operation;
    $dag->add_operation($build_adaptor_operation);
    $dag->connect_input(
        input_property => 'build_id',
        destination => $build_adaptor_operation,
        destination_property => 'build_id',
    );
    $dag->connect_input(
        input_property => 'variant_type',
        destination => $build_adaptor_operation,
        destination_property => 'variant_type',
    );
    return $build_adaptor_operation;
}

sub dag {
    #   Must return a Genome::WorkflowBuilder::DAG
    # these usually just consist of a build_adaptor
    # followed by a single command, but could be
    # more complex.

    # DAG INPUTS:
    #   build_id
    #   input_result  (Genome::SoftwareResult that has a
    #                  'output_file_path' accessor that refers
    #                  to a .vcf or .vcf.gz file. or a 'get_vcf'
    #                  accessor which takes a 'variant_type'
    #                  argument to refer to a .vcf or .vcf.gz
    #                  file)
    # DAG OUTPUTS:
    #   software_result (Same requirements as <input_result>)
    my $self = shift;
    return $self->_simple_dag;
}

sub _simple_dag {
    my $self = shift;

    my $dag = Genome::WorkflowBuilder::DAG->create(
        name => $self->name,
    );
    my $build_adaptor_operation = $self->connected_build_adaptor_operation($dag);

    my $run_operation = $self->connected_run_operation($dag);
    $self->_link(dag => $dag,
          adaptor => $build_adaptor_operation,
          target => $run_operation,
    );

    return $dag;
}

sub _link {
    my $self = shift;
    my %p = validate(@_, {
        dag => {isa => 'Genome::WorkflowBuilder::DAG'},
        adaptor => {isa => 'Genome::WorkflowBuilder::Command'},
        target => {isa => 'Genome::WorkflowBuilder::Command'},
    });

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
