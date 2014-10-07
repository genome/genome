package Genome::VariantReporting::Expert::Vep::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Expert::Vep::Run {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    has_input => [
        ensembl_version => {
            is => 'String',
        },
        custom_annotation_tags => {
            is => 'String',
            is_many => 1,
            is_optional => 1,
        },
        feature_list_ids => {
            is => 'HASH',
            doc => 'A hash keyed on INFO TAG with values of FeatureList IDs',
        },
        reference_fasta => {is => 'Path'},
        species => { is => 'Text', },
        plugins => {is => 'String',
                    is_many => 1,
                    is_optional => 1},
        joinx_version => { is => 'String',},
        plugins_version => {is => 'String',},
    ],
    has_param => [
        lsf_resource => {
            value => q{-R 'select[mem>32000] rusage[mem=32000]' -M 32000000},
        },
    ],
};

sub name {
    'vep';
}

sub result_class {
    'Genome::VariantReporting::Expert::Vep::RunResult';
}
