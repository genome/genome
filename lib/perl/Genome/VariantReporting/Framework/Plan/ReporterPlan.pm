package Genome::VariantReporting::Framework::Plan::ReporterPlan;

use strict;
use warnings FATAL => 'all';
use Genome;
use Set::Scalar;

class Genome::VariantReporting::Framework::Plan::ReporterPlan {
    is => 'Genome::VariantReporting::Framework::Plan::Base',
    has => [
        interpreter_plans => {
            is => 'Genome::VariantReporting::Framework::Plan::InterpreterPlan',
            is_many => 1,
        },
        filter_plans => {
            is => 'Genome::VariantReporting::Framework::Plan::FilterPlan',
            is_many => 1,
            is_optional => 1,
        },
    ],
};

sub category {
    'reporters';
}

sub children {
    my $self = shift;
    return ('interpreters' => [$self->interpreter_plans],
        'filters' => [$self->filter_plans]);
}

sub create_from_hashref {
    my $class = shift;
    my $name = shift;
    my $hashref = shift;

    my $self = $class->SUPER::create(
        name => $name,
        params => $hashref->{params},
    );

    my @filter_plans;
    for my $filter_name (keys %{$hashref->{filters}}) {
        push @filter_plans, Genome::VariantReporting::Framework::Plan::FilterPlan->create(
            name => $filter_name,
            params => $hashref->{filters}->{$filter_name},
        );
    }
    $self->filter_plans(\@filter_plans);

    my @interpreter_plans;
    for my $interpreter_name (keys %{$hashref->{interpreters}}) {
        push @interpreter_plans, Genome::VariantReporting::Framework::Plan::InterpreterPlan->create(
            name => $interpreter_name,
            params => $hashref->{interpreters}->{$interpreter_name},
        );
    }
    $self->interpreter_plans(\@interpreter_plans);

    return $self;
}

sub __plan_errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__plan_errors__;

    push @errors, $self->__interpreter_plan_errors__;

    return @errors;
}

sub __interpreter_plan_errors__ {
    my $self = shift;
    my $needed = Set::Scalar->new($self->get_class->requires_interpreters);
    my $have = Set::Scalar->new(map {$_->name} $self->interpreter_plans);

    my @errors;
    unless($needed->is_equal($have)) {
        if (my $still_needed = $needed - $have) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => [$still_needed->members],
                desc => sprintf("Interpreters required by reporter (%s) but not provided: (%s)", $self->name, join(",", $still_needed->members)),
            );
        }
        if (my $not_needed = $have - $needed) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => [$not_needed->members],
                desc => sprintf("Interpreters provided by plan but not required by reporter (%s): (%s)",
                    $self->name, join(",", $not_needed->members)),
            );
        }
    }
    return @errors;
}

sub requires_annotations {
    my $self = shift;
    my $needed = Set::Scalar->new();
    for my $plan ($self->interpreter_plans, $self->filter_plans) {
        $needed->insert($plan->get_class->requires_annotations);
    }
    return $needed->members;
}

1;
