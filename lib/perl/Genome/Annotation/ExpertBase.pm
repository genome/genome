package Genome::Annotation::ExpertBase;

use strict;
use warnings FATAL => 'all';
use Genome;
use Params::Validate qw(validate :types);

class Genome::Annotation::ExpertBase {
    is => 'Genome::Annotation::ComponentBase',
    is_abstract => 1,
};

sub adaptor_class {
    my $self = shift;
    my @parts = split(/::/, $self->class);
    pop @parts;
    my $result = join('::', @parts, 'Adaptor');
    if ($result->isa('Genome::Annotation::AdaptorBase')) {
        return $result;
    } else {
        die $self->error_message("Couldn't find an adaptor for expert (%s) at (%s)",
            $self->class, $result);
    }
}

sub build_adaptor_operation {
    my $self = shift;
    return Genome::WorkflowBuilder::Command->create(
        name => 'Get inputs from build',
        command => $self->adaptor_class,
    );
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
    die "Abstract";
}

sub _link {
    my $self = shift;
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
