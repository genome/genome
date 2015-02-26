package Genome::VariantReporting::Suite::Fpkm::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Fpkm::Run {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    has_planned_transient => [
        fpkm_file => {
            is => 'Path',
            is_translated => 1,
            doc => 'The fpkm file to use for annotation',
        },
        sample_name => {
            is => 'Text',
            is_translated => 1,
            doc => 'The sample to analyze',
        },
    ],
    doc => 'Annotate variants with FPKM value for the sample specified',
};

sub name {
    'fpkm';
}

sub result_class {
    'Genome::VariantReporting::Suite::Fpkm::RunResult';
}
