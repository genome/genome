package Genome::VariantReporting::Framework::Component::WithTranslatedInputs;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Framework::Utility;
use Exception::Class ('NoTranslationsException');

class Genome::VariantReporting::Framework::Component::WithTranslatedInputs {
    is => 'Genome::VariantReporting::Framework::Component::Base',
    attributes_have => [
        is_translated => {is => 'Boolean', default => 0},
    ],
};

sub translated_input_names {
    my $class = shift;
    return map {$_->property_name} $class->__meta__->properties(
        is_translated => 1,
        is_many => 0,
    );
}

sub translated_is_many_input_names {
    my $self = shift;
    return map {$_->property_name} $self->__meta__->properties(
        is_translated => 1,
        is_many => 1,
    );
}

1;
