package Genome::VariantReporting::Framework::Test::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Framework::Test::Adaptor {
    is => 'Genome::VariantReporting::Framework::Component::Adaptor',

    has_planned_output => [
        __planned__ => {},
    ],
    has_provided_output => [
        __provided__ => {},
    ],
};

sub name {
    "__test__";
}

1;
