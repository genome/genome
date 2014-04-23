package Genome::Annotation::Plan::Base;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Plan::Base {
    is => 'Genome::Annotation::ComponentBase',
    has => [
        name => {
            is => 'Text',
        },
        params => {
            is => 'HASH',
        },
    ],
};

sub category {
    die "Abstract (eg 'expert', 'filter', 'interpreter', or 'reporter')";
}

sub children {
    return ();
}

sub as_hashref {
    my $self = shift;

    my %body;
    if ($self->children) {
        $body{params} = $self->params;

        my %children = $self->children;
        for my $child_category (keys %children) {
            my @child_plans = @{$children{$child_category}};
            $body{$child_category} = {};
            for my $child_plan (@child_plans) {
                my $child_hashref = $child_plan->as_hashref;
                for my $key (keys %{$child_hashref}) {
                    $body{$child_category}{$key} = $child_hashref->{$key};
                }
            }
        }
    } else {
        for my $param_name (keys %{$self->params}) {
            $body{$param_name} = $self->params->{$param_name};
        }
    }

    my %result;
    $result{$self->name} = \%body;
    return \%result;
}

sub validate {
    my $self = shift;
    $self->validate_self;
    $self->validate_object;
}

sub validate_object {
    my $self = shift;
    if (my $object = $self->object) {
        $object->validate();
    }
}

sub validate_self {
    my $self = shift;
    $self->validate_params;
    $self->validate_children;
}

sub object {
    my $self = shift;
    return $self->factory->get_object($self->category,
        $self->name, $self->params);
}

sub factory {
    my $self = shift;
    return Genome::Annotation::Factory->create();
}

sub validate_children {
    my $self = shift;
    my %child_hash = $self->children;
    for my $children_of_category (values %child_hash) {
        for my $child (@{$children_of_category}) {
            $child->validate;
        }
    }
}

1;
