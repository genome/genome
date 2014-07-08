package Genome::VariantReporting::Vep::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Vep::Run {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    has_input => [
        ensembl_version => {
            is => 'String',
        },
        feature_list_ids_and_tags => {
            is => 'String',
            is_many => 1,
            doc => 'List of feature lists to be annotated in the 
                    INFO field, along with the tag to be used
                    e.g. 12345:SEGDUP,58676:ROI
                    The id and tag should be separated by a colon',
            is_optional => 1,
        },
        reference_build => {is => 'Genome::Model::Build::ReferenceSequence'},
        species => { is => 'Text', },
        terms => {is => 'String',},
        plugins => {is => 'String',
                    is_many => 1,
                    is_optional => 1},
        joinx_version => { is => 'String',},
        plugins_version => {is => 'String',},
    ],
    has_param => [
        lsf_resource => {
            value => q{-R 'select[mem>16000] rusage[mem=16000]' -M 16000000},
        },
    ],
};

sub name {
    'vep';
}

sub result_class {
    'Genome::VariantReporting::Vep::RunResult';
}
