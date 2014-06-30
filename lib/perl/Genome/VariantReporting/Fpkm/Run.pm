package Genome::VariantReporting::Fpkm::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Fpkm::Run {
    is => 'Genome::VariantReporting::Component::Expert::Command',
    has_input => [
        fpkm_file => {is => 'Path'},
        tumor_sample_name => {is => 'Text'},
    ],
};

sub name {
    'fpkm';
}

sub result_class {
    'Genome::VariantReporting::Fpkm::RunResult';
}
