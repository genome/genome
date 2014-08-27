package Genome::VariantReporting::FlankingRegions::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::FlankingRegions::Run {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    has_input => [
        flank_size => {
            is => 'String',
        },
        reference_fasta => {is => 'Path'},
        tumor_sample_name => {is => 'Text'},
    ],
};

sub name {
    'flanking-regions';
}

sub result_class {
    'Genome::VariantReporting::FlankingRegions::RunResult';
}
