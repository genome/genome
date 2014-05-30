package Genome::VariantReporting::Dag;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::WorkflowBuilder::DAG;
use Genome::WorkflowBuilder::Command;
use Genome::Utility::List qw(in);
use Params::Validate qw(validate validate_pos :types);

use Exporter 'import';

our @EXPORT_OK = qw(
    generate_dag
);

sub generate_dag {
    my ($plan, $variant_type) = validate_pos(@_, 1, 1);

    my $dag = Genome::WorkflowBuilder::DAG->create(
        name => sprintf('Variant Reporting (%s)', $variant_type),
    );

    my $last_expert_op = connect_experts($dag, $plan);
    connect_report_generator($dag, $last_expert_op);

    return $dag;
}

sub connected_initial_operation {
    my $dag = shift;

    my $op = Genome::WorkflowBuilder::Command->create(
        name => 'Get Initial Vcf Result from Build',
        command => 'Genome::VariantReporting::GetInitialVcfResult',
    );
    connect_to_dag(
        dag => $dag,
        target => $op,
    );
    return $op;
}

sub connect_experts {
    my $dag = shift;
    my $plan = shift;

    my $last_op = connected_initial_operation($dag);

    for my $expert (experts($plan)) {
        my $expert_op = $expert->dag;
        connect_to_previous(
            dag => $dag,
            previous => $last_op,
            target => $expert_op,
        );
        connect_to_dag(
            dag => $dag,
            target => $expert_op,
        );
        $last_op = $expert_op;
    }
    return $last_op;
}

sub connect_to_dag {
    my %p = validate(@_, {
        dag => {isa => 'Genome::WorkflowBuilder::DAG'},
        target => {type => OBJECT},
    });

    $p{dag}->add_operation($p{target});
    for my $name qw(build_id variant_type plan_json) {
        if (in($name, $p{target}->input_properties)) {
            $p{dag}->connect_input(
                input_property => $name,
                destination => $p{target},
                destination_property => $name,
            );
        }
    }
}

sub connect_to_previous {
    my %p = validate(@_, {
        dag => {isa => 'Genome::WorkflowBuilder::DAG'},
        previous => {type => OBJECT | UNDEF},
        target => {type => OBJECT},
    });

    $p{dag}->create_link(
        source => $p{previous},
        source_property => 'output_result',
        destination => $p{target},
        destination_property => 'input_result',
    );
}

sub experts {
    my $plan = shift;

    my @unsorted_experts = map {$_->object} $plan->expert_plans;

    # sort first by higher priority then by alphabetically earlier name.
    return sort {
        $b->priority <=> $a->priority ||
        $a->name cmp $b->name
    } @unsorted_experts;
}

sub connect_report_generator {
    my ($dag, $last_expert_op) = validate_pos(@_, 1, 1);

    my $report_generator_op = Genome::WorkflowBuilder::Command->create(
        name => 'Generate Reports',
        command => 'Genome::VariantReporting::ReportGeneratorWrapper',
    );
    connect_to_previous(
        dag => $dag,
        previous => $last_expert_op,
        target => $report_generator_op,
    );
    connect_to_dag(
        dag => $dag,
        target => $report_generator_op,
    );

    $dag->connect_input(
        input_property => 'output_directory',
        destination => $report_generator_op,
        destination_property => 'output_directory',
    );
    $dag->connect_output(
        output_property => 'output_directory',
        source => $report_generator_op,
        source_property => 'output_directory',
    );
}


1;
