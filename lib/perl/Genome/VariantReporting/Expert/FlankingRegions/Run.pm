package Genome::VariantReporting::Expert::FlankingRegions::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Expert::FlankingRegions::Run {
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
    'Genome::VariantReporting::Expert::FlankingRegions::RunResult';
}
