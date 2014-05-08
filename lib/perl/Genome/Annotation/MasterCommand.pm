package Genome::Annotation::MasterCommand;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::WorkflowBuilder::DAG;
use Genome::Utility::List qw(in);
use Genome::WorkflowBuilder::Command;
use Params::Validate qw(validate validate_pos :types);

class Genome::Annotation::MasterCommand {
    is => 'Command::V2',
    has_input => [
        build => {
            is => 'Genome::Model::Build',
            is_optional => 1,
        },
        variant_type => {
            is => 'Text',
            is_output => 1,
            valid_values => ['snvs', 'indels'],
        },
        output_directory => {
            is => 'Path',
            is_optional => 1,
            is_output => 1,
        },
        log_directory => {
            is => 'Path',
            is_optional => 1,
        },
        plan => {
            is => 'Genome::Annotation::Plan',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->status_message("Constructing workflow from plan.");
    my $dag = $self->dag;

    Genome::Sys->create_directory($self->log_directory);
    $dag->log_dir($self->log_directory);

    $self->status_message("Executing workflow.");
    $dag->execute(
        build_id => $self->build->id,
        variant_type => $self->variant_type,
        output_directory => $self->output_directory,
        plan_json => $self->plan->as_json,
    );
    return 1;
}

sub dag {
    my $self = shift;

    my $dag = Genome::WorkflowBuilder::DAG->create(
        name => sprintf('Annotation (%s)', $self->variant_type),
    );

    my $last_expert_op = $self->connect_experts($dag);
    $self->connect_report_generator($dag, $last_expert_op);

    return $dag;
}

sub connected_initial_operation {
    my ($self, $dag) = @_;

    my $op = Genome::WorkflowBuilder::Command->create(
        name => 'Get Initial Vcf Result from Build',
        command => 'Genome::Annotation::GetInitialVcfResult',
    );
    $self->connect_to_dag(
        dag => $dag,
        target => $op,
    );
    return $op;

}

sub connect_experts {
    my $self = shift;
    my $dag = shift;

    my $last_op = $self->connected_initial_operation($dag);

    for my $expert ($self->experts) {
        my $expert_op = $expert->dag;
        $self->connect_to_previous(
            dag => $dag,
            previous => $last_op,
            target => $expert_op,
        );
        $self->connect_to_dag(
            dag => $dag,
            target => $expert_op,
        );
        $last_op = $expert_op;
    }
    return $last_op;
}

sub connect_to_dag {
    my $self = shift;
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
    my $self = shift;
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
    my $self = shift;

    my @unsorted_experts = map {$_->object} $self->plan->expert_plans;

    # sort first by higher priority then by alphabetically earlier name.
    return sort {
        $b->priority <=> $a->priority ||
        $a->name cmp $b->name
    } @unsorted_experts;
}

sub connect_report_generator {
    my ($self, $dag, $last_expert_op) = validate_pos(@_, 1, 1, 1);

    my $report_generator_op = Genome::WorkflowBuilder::Command->create(
        name => 'Generate Reports',
        command => 'Genome::Annotation::ReportGeneratorWrapper',
    );
    $self->connect_to_previous(
        dag => $dag,
        previous => $last_expert_op,
        target => $report_generator_op,
    );
    $self->connect_to_dag(
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
