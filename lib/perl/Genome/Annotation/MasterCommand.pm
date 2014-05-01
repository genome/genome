package Genome::Annotation::MasterCommand;

use strict;
use warnings FATAL => 'all';
use Genome;
use Params::Validate qw(validate validate_pos :types);

class Genome::Annotation::MasterCommand {
    is => 'Command::V2',
    has_input => [
        build_id => {
            is => 'Text',
        },
        variant_type => {
            is => 'Text',
            is_output => 1,
            valid_values => ['snvs', 'indels'],
        },
        output_directory => {
            is => 'Path',
            is_output => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    $self->dag->execute(
        build_id => $self->build_id,
        variant_type => $self->variant_type,
    );
    return 1;
}

sub dag {
    my $self = shift;

    my $dag = Genome::WorkflowBuilder::DAG->create(
        name => 'Annotation',
    );

    my $last_expert_op = $self->connect_experts($dag);
    $self->connect_report_generator($dag, $last_expert_op);

    return $dag;
}

sub connect_experts {
    my $self = shift;
    my $dag = shift;

    my $last_expert_op;
    for my $expert ($self->experts) {
        my $expert_op = $expert->dag;
        $self->connect_to_previous(
            dag => $dag,
            previous => $last_expert_op,
            target => $expert_op,
        );
        $last_expert_op = $expert_op;
    }
    return $last_expert_op;
}

sub connect_to_previous {
    my $self = shift;
    my %p = validate(@_, {
        dag => {isa => 'Genome::WorkflowBuilder::DAG'},
        previous => {type => OBJECT | UNDEF},
        target => {type => OBJECT},
    });

    $p{dag}->add_operation($p{target});

    if (defined $p{previous}) {
        $p{dag}->create_link(
            source => $p{previous},
            source_property => 'output_result',
            destination => $p{target},
            destination_property => 'input_result',
        );
    } else {
        $p{dag}->connect_input(
            input_property => 'initial_vcf_result',
            destination => $p{target},
            destination_property => 'input_result',
        );
    }

    for my $name qw(build_id variant_type) {
        $p{dag}->connect_input(
            input_property => $name,
            destination => $p{target},
            destination_property => $name,
        );
    }
}

sub build {
    my $self = shift;

    my $build = Genome::Model::Build->get($self->build_id);
    if ($build) {
        return $build;
    } else {
        die $self->error_message("Couldn't find a build for id (%s)",
            $self->build_id);
    }
}

sub plan {
    my $self = shift;

    return $self->build->annotation_plan($self->variant_type);
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
