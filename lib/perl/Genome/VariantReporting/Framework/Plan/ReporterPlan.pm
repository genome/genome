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
    while (my ($name, $params) = each %{$hashref->{filters}}) {
        push @filter_plans, Genome::VariantReporting::Framework::Plan::FilterPlan->create(
            name => $name,
            params => $params,
        );
    }
    $self->filter_plans(\@filter_plans);

    my @interpreter_plans;
    while (my ($name, $params) = each %{$hashref->{interpreters}}) {
        push @interpreter_plans, Genome::VariantReporting::Framework::Plan::InterpreterPlan->create(
            name => $name,
            params => $params,
        );
    }
    $self->interpreter_plans(\@interpreter_plans);

    return $self;
}

sub __translation_errors__ {
    my ($self, $provider) = @_;
    my @errors;
    push @errors, $self->object->_translation_errors($provider->translations, $self->object->name);
    for my $child (values(%{$self->object->interpreters}), values(%{$self->object->filters})) {
        push @errors, $child->_translation_errors($provider->translations, $child->name);
    }
    return @errors;
}

sub __plan_errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__plan_errors__;

    push @errors, $self->__interpreter_plan_errors__;
    push @errors, $self->__filter_plan_errors__;

    return @errors;
}

sub __interpreter_plan_errors__ {
    my $self = shift;
    my $needed = Set::Scalar->new($self->get_class->requires_interpreters);
    my $have = Set::Scalar->new(map {$_->name} $self->interpreter_plans);

    my @errors;
    if (my $still_needed = $needed - $have) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => [$still_needed->members],
            desc => sprintf("Interpreters required by reporter (%s) but not provided: (%s)", $self->name, join(",", $still_needed->members)),
        );
    }
    return @errors;
}

sub __filter_plan_errors__ {
    my $self = shift;

    unless ($self->get_class->allows_hard_filters) {
        if ($self->filter_plans) {
            return UR::Object::Tag->create(
                type => 'error',
                desc => sprintf("Reporter (%s) does not allow any hard filters.  Move them to the interpreters section if they should be soft filters. (%s)",
                    $self->name, join(",", map {$_->name} $self->filter_plans)),
            );
        }
    }
    return;
}

sub requires_annotations {
    my $self = shift;
    my $needed = Set::Scalar->new();
    for my $plan ($self->interpreter_plans, $self->filter_plans) {
        $needed->insert($plan->get_class->requires_annotations);
    }
    return $needed->members;
}

sub object {
    my $self = shift;

    my $reporter = $self->SUPER::object(@_);

    for my $filter_plan ($self->filter_plans) {
        $reporter->add_filter_object($filter_plan->object(@_));
    }
    for my $interpreter_plan ($self->interpreter_plans) {
        $reporter->add_interpreter_object($interpreter_plan->object(@_));
    }

    return $reporter;
}
Memoize::memoize("object", LIST_CACHE => 'MERGE');


1;
