package Genome::VariantReporting::Framework::Component::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;
use Params::Validate qw(validate_pos validate :types);
use Genome::WorkflowBuilder::DAG;
use Genome::WorkflowBuilder::Command;
use Carp qw(confess);

class Genome::VariantReporting::Framework::Component::Expert {
    is => 'Genome::VariantReporting::Framework::Component::Base',
    is_abstract => 1,
};

sub name {
    my $self = shift;
    my $class = $self->class;
    confess "Abstract method 'name' must be defined in class '$class'";
}

sub provides_annotations {
    my $class = shift;
    return ('variant_calls', $class->name);
}

sub priority {
    # Higher priority experts are run earlier.
    # Ties are broken by expert name (alphabetically early names are
    # run earlier).
    return 0;
}

sub run_class {
    my $self = shift;
    my $factory = Genome::VariantReporting::Framework::Factory->create();
    return $factory->get_class('runners', $self->name);
}

sub adaptor_class {
    my $self = shift;
    my $factory = Genome::VariantReporting::Framework::Factory->create();
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
        input_property => 'input_vcf',
        destination => $run_operation,
        destination_property => 'input_vcf',
    );
    for my $name qw(output_result output_vcf) {
        $dag->connect_output(
            output_property => $name,
            source => $run_operation,
            source_property => $name,
        );
    }
    return $run_operation;
}

sub adaptor_operation {
    my $self = shift;
    return Genome::WorkflowBuilder::Command->create(
        name => 'Get inputs from provider and plan',
        command => $self->adaptor_class,
    );
}

sub connected_adaptor_operation {
    my ($self, $dag) = validate_pos(@_, 1, 1);

    my $adaptor_operation = $self->adaptor_operation;
    $dag->add_operation($adaptor_operation);
    for my $name qw(provider_json variant_type plan_json) {
        $dag->connect_input(
            input_property => $name,
            destination => $adaptor_operation,
            destination_property => $name,
        );
    }
    return $adaptor_operation;
}

sub dag {
    #   Must return a Genome::WorkflowBuilder::DAG
    # these usually just consist of an adaptor
    # followed by a single command, but could be
    # more complex.

    # DAG INPUTS:
    #   input_vcf
    #   provider_json
    #   variant_type
    #   plan_json
    #
    # DAG OUTPUTS:
    #   output_vcf
    my $self = shift;
    return $self->_simple_dag;
}

sub _simple_dag {
    my $self = shift;

    my $dag = Genome::WorkflowBuilder::DAG->create(
        name => $self->name,
    );
    my $adaptor_operation = $self->connected_adaptor_operation($dag);

    my $run_operation = $self->connected_run_operation($dag);
    $self->_link(dag => $dag,
          adaptor => $adaptor_operation,
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
        next if $name eq 'input_vcf';
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
