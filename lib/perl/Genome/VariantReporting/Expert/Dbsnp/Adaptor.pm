package Genome::VariantReporting::Expert::Dbsnp::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Expert::Dbsnp::Adaptor {
    is => 'Genome::VariantReporting::AdaptorBase',

    has_planned_output => [
        joinx_version => { is  => 'Version', },
        info_string => { is => 'Text', },
    ],
    has_output => [
        known_variants => {
            is => 'Genome::Model::Build::ImportedVariationList',
        },
    ],
};

sub name {
    "dbsnp";
}

sub resolve_expert_specific_attributes_from_build {
    my $self = shift;

    $self->known_variants($self->build->previously_discovered_variations_build);
    return;
}

1;
