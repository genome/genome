package Genome::VariantReporting::Framework::Test::WithTranslationsExpert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Framework::Test::WithTranslationsExpert {
    is => 'Genome::VariantReporting::Framework::Component::Expert',
};

sub name {
    '__translated_test__';
}

1;
