package Genome::Annotation::Expert::Dbsnp::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Expert::Dbsnp::Adaptor {
    is => 'Genome::Annotation::AdaptorBase',

    has_planned_output => [
        joinx_version => { is  => 'Version', },
        info_string => { is => 'Text', },
    ],
    has_output => [
        known_variants => {
            is => 'Genome::Model::Build::ImportedVariationList',
            is_many => 1,
        },
    ],
};

sub name {
    "dbsnp";
}

sub resolve_expert_specific_attributes_from_build {
    my $self = shift;

    $self->known_variants([$self->build->previously_discovered_variations_build,]);
    return;
}

1;
