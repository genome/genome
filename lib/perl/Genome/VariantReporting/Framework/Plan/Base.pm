package Genome::VariantReporting::Framework::Plan::Base;

use strict;
use warnings FATAL => 'all';
use Genome;
use Memoize qw();
use Params::Validate qw(validate_pos :types);
use Data::Dump qw(pp);
use Exception::Class ('NoTranslationsException');

class Genome::VariantReporting::Framework::Plan::Base {
    is => 'Genome::VariantReporting::Framework::Component::Base',
    has => [
        name => {
            is => 'Text',
        },
        params => {
            is => 'HASH',
            default => {},
        },
    ],
};

our $FACTORY = Genome::VariantReporting::Framework::Factory->create();

sub category {
    die "Abstract (eg 'expert', 'filter', 'interpreter', or 'report')";
}

sub children {
    return ();
}

sub as_hashref {
    my $self = shift;

    my %body;
    if ($self->children) {
        $body{params} = $self->params if keys %{$self->params};

        my %children = $self->children;
        while (my ($child_category, $child_plans_ref) = each %children) {
            my @child_plans = @$child_plans_ref;
            $body{$child_category} = {};
            for my $child_plan (@child_plans) {
                my $child_hashref = $child_plan->as_hashref;
                while (my ($name, $child) = each %$child_hashref) {
                    $body{$child_category}{$name} = $child;
                }
            }
        }
    } else {
        while (my ($name, $value) = each %{$self->params}) {
            $body{$name} = $value;
        }
    }

    my %result;
    $result{$self->name} = \%body;
    return \%result;
}

sub get_class {
    my $self = shift;
    return $FACTORY->get_class($self->category, $self->name);
}

sub object {
    my ($self) = validate_pos(@_, {type => OBJECT});

    return $FACTORY->get_object($self->category,
            $self->name, $self->params);
}
Memoize::memoize("object", LIST_CACHE => 'MERGE');

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

sub __translation_errors__ {
    my $self = shift;
    my $translations = shift;

    return Genome::VariantReporting::Framework::Utility::get_missing_errors($self->get_class->name,
        $translations, $self->needed_translations, "Translations", "component");
}

sub needed_translations {
    my $self = shift;
    return Set::Scalar->new(map {$self->params->{$_}} $self->get_class->required_translated_input_names);
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

sub translate {
    my ($self, $translations) = @_;
    return $self->_translate($translations, 'params', $self->get_class);
}

sub _translate {
    my ($self, $translations, $params_accessor, $object_class) = @_;

    for my $name ($object_class->all_translated_input_names) {
        my $old_value = $self->$params_accessor->{$name};
        next unless defined($old_value);
        $self->$params_accessor->{$name} = $self->_translate_single($old_value, $translations, $name);
    }

    for my $name ($object_class->all_translated_is_many_input_names) {
        my $planned_value = $self->$params_accessor->{$name};
        next unless defined($planned_value);

        my @old_values = @{$planned_value};
        my $data_type = $object_class->__meta__->property_meta_for_name($name)->data_type;
        if (defined($data_type) && $data_type eq 'ARRAY') {
            $self->$params_accessor->{$name} = [map {$self->_translate_single($_, $translations, $name)} @old_values];
        }
        else {
            my @translated;
            for my $value (@old_values) {
                my $translated_value = $self->_translate_single($value, $translations, $name);
                if (ref($translated_value) eq 'ARRAY') {
                    push @translated, @$translated_value;
                }
                else {
                    push @translated, $translated_value;
                }
            }
            my $is_optional = $object_class->__meta__->property_meta_for_name($name)->is_optional;
            if (scalar(@translated) or $is_optional) {
                $self->$params_accessor->{$name} = \@translated;
            } else {
                die sprintf("No translations for input (%s) for plan (%s) which is a (%s):\n%s\nTranslations:\n%s\n",
                    $name, $self->name, $self->category, pp($self->as_hashref), pp($translations));

            }
        }
    }
    return;
}

sub _translate_single {
    my ($self, $old_value, $translations, $name) = @_;

    unless (defined($translations)) {
        NoTranslationsException->throw(
            error => sprintf(
                "Could not translate input (%s) with value (%s) for plan (%s) which is a (%s). No translations provided.",
                $name, $old_value, $self->name, $self->category,
            )
        );
    }

    unless (defined($old_value)) {
        die sprintf("No old_value for input (%s) for plan (%s) which is a (%s):\n%s",
            $name, $self->name, $self->category, pp($self->as_hashref));
    }
    if (exists($translations->{$old_value})) {
        return $translations->{$old_value};
    } else {
        die sprintf("Could not translate input (%s) with value (%s) for plan (%s) which is a (%s). Available translations are: %s",
            $name, $old_value, $self->name, $self->category, Data::Dumper::Dumper($translations));
    }
}


1;
