package Genome::VariantReporting::Suite::FlankingRegions::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::FlankingRegions::Adaptor {
    is => "Genome::VariantReporting::Framework::Component::Adaptor",

    has_planned_output => [
        flank_size => {
            is => 'Integer',
        },
    ],
    has_provided_output => [
        reference_fasta => {
            is => 'Path',
        },
        tumor_sample_name => {
            is => 'Text',
        },
    ],
};

sub name {
    'flanking-regions';
}

1;
