package Genome::VariantReporting::Vep::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Vep::Adaptor {
    is => "Genome::VariantReporting::Framework::Component::Adaptor",

    has_planned_output => [
        ensembl_version => {
            is => 'String',
        },
        species => { is => 'Text', },
        terms => {is => 'String',},
        plugins => {is => 'String',
                    is_many => 1,
                    is_optional => 1},
        custom_annotation_tags => {
            is => 'String',
            is_many => 1,
            is_optional => 1,
        },
        joinx_version => { is => 'String', },
        plugins_version => { is => 'String', },
    ],
    has_provided_output => [
        feature_list_ids => {
            is => 'HASH',
            doc => 'A hash keyed on INFO TAG with values of FeatureList IDs',
        },
        reference_fasta => {
            is => 'Path',
        },
    ],
};

sub name {
    'vep';
}

1;
