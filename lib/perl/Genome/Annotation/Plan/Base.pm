package Genome::Annotation::Plan::Base;

use strict;
use warnings FATAL => 'all';
use Genome;
use Memoize qw(memoize);

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

# We overide ComponentBase validate because only plans have objects to validate
sub validate {
    my $self = shift;

    $self->SUPER::validate;
    $self->validate_object;
}

sub validate_object {
    my $self = shift;
    if (my $object = $self->object) {
        $object->validate();
    }
}

sub get_class {
    my $self = shift;
    return $self->factory->get_class($self->category, $self->name);
}

sub object {
    my $self = shift;
    my %overrides = @_;
    return $self->factory->get_object($self->category,
            $self->name, $self->params, \%overrides);
}
Memoize::memoize("object");

sub factory {
    my $self = shift;
    return Genome::Annotation::Factory->create();
}
Memoize::memoize("factory");

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;

    my @child_errors = $self->__child_errors__;
    push @errors, @child_errors;

    my @class_errors = $self->__class_errors__;
    push @errors, @class_errors;

    # We can't check plan or object errors if there are class problems
    # i.e. the class does not exist
    if (@child_errors or @class_errors) {
        return @errors;
    } else {
        push @errors, $self->__object_errors__, $self->__plan_errors__;
    }

    return @errors;
}

sub __child_errors__ {
    my $self = shift;

    my @errors;
    my %child_hash = $self->children;
    for my $children_of_category (values %child_hash) {
        for my $child (@{$children_of_category}) {
            push @errors, $child->__errors__;
        }
    }

    return @errors;
}

sub __plan_errors__ {
    return;
}

sub __object_errors__ {
    my $self = shift;

    my $object;
    eval {
        $object = $self->object;
    };

    if (my $error = $@) {
        return UR::Object::Tag->create(
            type => 'error',
            properties => [],
            desc => $error,
        );
    } else {
        return $object->__errors__;
    }
}

sub __class_errors__ {
    my $self = shift;
    my @errors;
    eval {
        my $class = $self->get_class;
    };

    if (my $error = $@) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => [],
            desc => $error,
        );
    }

    return @errors;
}

1;
