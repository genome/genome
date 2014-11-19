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
        reference_fasta => {
            is => 'Path',
            is_translated => 1,
        },
        tumor_sample_name => {
            is => 'Text',
            is_translated => 1,
        },
    ],
    doc => 'Extract flanking sequence around variant and reference allele',
};

sub name {
    'flanking-regions';
}

1;
