package Genome::VariantReporting::Suite::BamReadcount::Annotate;

use strict;
use warnings FATAL => 'all';

class Genome::VariantReporting::Suite::BamReadcount::Annotate {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    has_input => [
        readcount_results => {
            is => 'Genome::VariantReporting::Suite::BamReadcount::RunResult',
            is_many => 1,
        },
    ],
    has_optional_output => [
        output_vcf => {
            is => 'Path',
        },
    ],
};

sub result_class {
    'Genome::VariantReporting::Suite::BamReadcount::AnnotateResult';
}

sub resolve_plan_attributes {
    return;
}

1;
