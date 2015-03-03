package Genome::VariantReporting::Suite::BamReadcount::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::WorkflowBuilder::DAG;
use Genome::WorkflowBuilder::Command;
use Params::Validate qw(validate_pos);

class Genome::VariantReporting::Suite::BamReadcount::Expert {
    is => 'Genome::VariantReporting::Framework::Component::Expert',
};

sub name {
    'bam-readcount';
}

sub dag {
    my $self = shift;

    my $dag = Genome::WorkflowBuilder::DAG->create(
        name => 'BamReadcount',
    );
    my $adaptor_op = $self->connected_adaptor_operation($dag);

    my $run_op = $self->run_op;
    $dag->add_operation($run_op);
    $run_op->parallel_by('aligned_bam_result_id');
    for my $name (qw(process_id input_vcf variant_type plan_json)) {
        $dag->connect_input(
            input_property => $name,
            destination => $run_op,
            destination_property => $name,
        );
    }
    $self->_link(dag => $dag,
          adaptor => $adaptor_op,
          target => $run_op,
    );

    my $annotate_op = $self->annotate_op;
    $dag->add_operation($annotate_op);
    for my $name (qw(process_id input_vcf variant_type plan_json)) {
        $dag->connect_input(
            input_property => $name,
            destination => $annotate_op,
            destination_property => $name,
        );
    }
    $dag->create_link(
        source => $run_op,
        source_property => 'output_result',
        destination => $annotate_op,
        destination_property => 'readcount_results',
    );

    $dag->connect_output(
        output_property => 'output_vcf',
        source => $annotate_op,
        source_property => 'output_vcf',
    );

    return $dag;
}

sub run_op {
    return Genome::WorkflowBuilder::Command->create(
        name => 'Run bam-readcount',
        command => 'Genome::VariantReporting::Suite::BamReadcount::Run',
    );
}

sub annotate_op {
    return Genome::WorkflowBuilder::Command->create(
        name => 'Annotate vcf with readcounts',
        command => 'Genome::VariantReporting::Suite::BamReadcount::Annotate',
    );
}

sub adaptor_operation {
    my $self = shift;
    return Genome::WorkflowBuilder::Command->create(
        name => 'Get inputs from plan',
        command => 'Genome::VariantReporting::Suite::BamReadcount::Adaptor',
    );
}

sub connected_adaptor_operation {
    my ($self, $dag) = validate_pos(@_, 1, 1);

    my $adaptor_operation = $self->adaptor_operation;
    $dag->add_operation($adaptor_operation);
    $dag->connect_input(
        input_property => 'plan_json',
        destination => $adaptor_operation,
        destination_property => 'plan_json',
    );
    return $adaptor_operation;
}

sub run_class {
    return 'Genome::VariantReporting::Suite::BamReadcount::Run';
}

1;
