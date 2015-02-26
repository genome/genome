package Genome::VariantReporting::Framework::Test::WithTranslationsRun;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Framework::Test::WithTranslationsRun {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    has_planned_transient => [
        __planned__ => {
            is_translated => 1,
        },
        __input__ => {
            is_many => 1,
            is_translated => 1,
        },
    ],
};

sub name {
    '__translated_test__';
}

sub result_class {
    'Genome::VariantReporting::Framework::Test::WithTranslationsRunResult';
}
