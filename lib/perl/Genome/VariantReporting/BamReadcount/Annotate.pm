package Genome::VariantReporting::BamReadcount::Annotate;

use strict;
use warnings FATAL => 'all';

class Genome::VariantReporting::BamReadcount::Annotate {
    is => 'Genome::VariantReporting::Component::Expert::Command',
    has_input => [
        readcount_results => {
            is => 'Genome::VariantReporting::BamReadcount::RunResult',
            is_many => 1,
        },
    ],
};

sub result_class {
    'Genome::VariantReporting::BamReadcount::AnnotateResult';
}

1;
