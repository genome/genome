package Genome::VariantReporting::Framework::Plan::Base;

use strict;
use warnings FATAL => 'all';
use Genome;
use Memoize qw(memoize);

class Genome::VariantReporting::Framework::Plan::Base {
    is => 'Genome::VariantReporting::Framework::Component::Base',
    has => [
        name => {
            is => 'Text',
        },
        params => {
            is => 'HASH',
        },
    ],
};

our $FACTORY = Genome::VariantReporting::Framework::Factory->create();

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
        # If we don't have any params, keys will crash
        if (defined $self->params) {
            for my $param_name (keys %{$self->params}) {
                $body{$param_name} = $self->params->{$param_name};
            }
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
    return $FACTORY->get_class($self->category, $self->name);
}

# We want to be able to get different perl objects when object vs
# object_with_translations is called -> don't memoize this
sub _object {
    my $self = shift;
    my %overrides = @_;
    return $FACTORY->get_object($self->category,
            $self->name, $self->params, \%overrides);
}

sub object {
    my $self = shift;
    my %overrides = @_;
    return $self->_object(%overrides);
}
Memoize::memoize("object", LIST_CACHE => 'MERGE');

sub object_with_translations {
    my $self = shift;
    my $translations = shift;
    my %overrides = @_;

    # Get a new perl object
    # For child classes that overwrite object: Using $self->_object rather
    # than $self->object ensures that we always call the base class' function
    my $object = $self->_object();

    if ($object->can('translate_inputs')) {
        $object->translate_inputs($translations);
    }
    return $object;
}
Memoize::memoize("object_with_translations", LIST_CACHE => 'MERGE');

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
            desc => sprintf("Problems with the plan for name (%s) of category (%s)", $self->name, $self->category) . "\n$error",
        );
    } else {
        my @errors = $object->__errors__;
        if (@errors) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                desc => sprintf("Problems with the plan for name (%s) of category (%s)", $self->name, $self->category),
            );
        }
        return @errors;
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
