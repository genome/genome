package Genome::VariantReporting::Framework::Component::Reporter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Carp qw(confess);

class Genome::VariantReporting::Framework::Component::Reporter {
    is => [
        'Genome::VariantReporting::Framework::Component::Base',
        'Genome::VariantReporting::Framework::Component::WithTranslatedInputs',
    ],
    is_abstract => 1,
    has_transient_optional => [
        filters => {
            is => 'HASH',
        },
        interpreters => {
            is => 'HASH',
        },
    ],
};

sub name {
    my $self = shift;
    my $class = $self->class;
    confess "Abstract method 'name' must be defined in class '$class'";
}

sub requires_interpreters {
    my $self = shift;
    my $class = $self->class;
    confess "Abstract method 'requires_interpreters' must be defined in class '$class'";
}

sub allows_hard_filters {
    return 1;
}

sub initialize {
    # this gets called before interpretations are given to ->report method
    return;
}

sub report {
    my ($self, $interpretations) = shift;
    # do something with the interpretations to produce one or more reports
    return;
}

sub finalize {
    # this gets called after interpretations are given to ->report method
    return;
}

sub add_filter_object {
    my ($self, $filter) = @_;

    my $filters_ref = $self->filters || {};
    $filters_ref->{$filter->name} = $filter;
    $self->filters($filters_ref);
}

sub add_interpreter_object {
    my ($self, $interpreter) = @_;

    my $interpreters_ref = $self->interpreters || {};
    $interpreters_ref->{$interpreter->name} = $interpreter;
    $self->interpreters($interpreters_ref);
}

1;
