package Genome::VariantReporting::Framework::Test::WithTranslationsAdaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Framework::Test::WithTranslationsAdaptor {
    is => 'Genome::VariantReporting::Framework::Component::Adaptor',

    has_planned_translated_output => [
        __planned__ => {},
        __provided__ => {},
    ],
};

sub name {
    "__translated_test__";
}

1;
