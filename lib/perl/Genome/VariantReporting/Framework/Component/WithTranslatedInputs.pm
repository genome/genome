package Genome::VariantReporting::Framework::Component::WithTranslatedInputs;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Framework::Utility;
use Exception::Class ('NoTranslationsException');

class Genome::VariantReporting::Framework::Component::WithTranslatedInputs {
    attributes_have => [
        is_translated => {is => 'Boolean', default => 0},
    ],
};

sub translate_inputs {
    my ($self, $translations) = @_;
    for my $name ($self->translated_input_names) {
        my $old_value = $self->$name;
        $self->$name($self->translate($old_value, $translations, $name));
    }

    for my $name ($self->translated_is_many_input_names) {
        my @old_values = $self->$name;
        if ($self->__meta__->property_meta_for_name($name)->data_type eq 'ARRAY') {
            $self->$name([map {$self->translate($_, $translations, $name)} @old_values]);
        }
        else {
            my @translated;
            for my $value (@old_values) {
                my $translated_value = $self->translate($value, $translations, $name);
                if (ref($translated_value) eq 'ARRAY') {
                    push @translated, @$translated_value;
                }
                else {
                    push @translated, $translated_value;
                }
            }
            $self->$name(\@translated);
        }
    }
}

sub translated_input_names {
    my $self = shift;
    return map {$_->property_name} $self->__meta__->properties(
        is_translated => 1,
        is_many => 0,
    );
}

sub needed_translations {
    my $self = shift;
    return Set::Scalar->new(map {$self->$_} $self->translated_input_names);
}

sub _translation_errors {
    my ($self, $translations, $component_name) = @_;
    return Genome::VariantReporting::Framework::Utility::get_missing_errors($component_name,
        $translations, $self->needed_translations, "Translations", "component");
}

sub translated_is_many_input_names {
    my $self = shift;
    return map {$_->property_name} $self->__meta__->properties(
        is_translated => 1,
        is_many => 1,
    );
}

sub translate {
    my ($self, $old_value, $translations, $name) = @_;

    unless (defined($translations)) {
        NoTranslationsException->throw(
            error => sprintf(
                "Could not translate input (%s) with value (%s) for object (%s). No translations provided.",
                $name, $old_value, $self->class
            )
        );
    }

    if (exists($translations->{$old_value})) {
        return $translations->{$old_value};
    } else {
        die sprintf("Could not translate input (%s) with value (%s) for object (%s). Available translations are: %s",
            $name, $old_value, $self->class, Data::Dumper::Dumper($translations));
    }
}

1;
