package Genome::VariantReporting::Suite::BamReadcount::Annotate;

use strict;
use warnings FATAL => 'all';

class Genome::VariantReporting::Suite::BamReadcount::Annotate {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    attributes_have => {
        is_argument => {
            is => "Boolean",
            default => 0,
        },
    },
    has_input => [
        input_vcf => {
            is => 'Path',
            is_argument => 1,
        },
        variant_type => {
            is => 'Text',
            is_argument => 1,
            valid_values => ['snvs', 'indels'],
            doc => "The type of variant the input_result represents",
        },
        process_id => {
            is => 'Text',
        },
    ],
    has_input_argument => [
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
    has_transient_optional => [
        requestor => {
            is => 'Genome::Process',
            id_by => 'process_id',
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
