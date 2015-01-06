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

sub all_translated_inputs {
    my $class = shift;
    return $class->__meta__->properties(
        is_translated => 1,
        is_many => 0,
    );
}

sub all_translated_is_many_inputs {
    my $class = shift;
    return $class->__meta__->properties(
        is_translated => 1,
        is_many => 1,
    );
}

sub all_translated_input_names {
    my $class = shift;
    return map {$_->property_name} $class->all_translated_inputs;
}

sub all_translated_is_many_input_names {
    my $class = shift;
    return map {$_->property_name} $class->all_translated_is_many_inputs;
}

sub required_translated_input_names {
    my $class = shift;
    return map {$_->property_name} grep {!$_->is_optional}
        $class->all_translated_inputs;
}

sub required_translated_is_many_input_names {
    my $class = shift;
    return map {$_->property_name} grep {!$_->is_optional}
        $class->all_translated_is_many_inputs;
}

1;
