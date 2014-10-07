package Genome::VariantReporting::Expert::Fpkm::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;
use File::Spec;

class Genome::VariantReporting::Expert::Fpkm::Adaptor {
    is => "Genome::VariantReporting::Framework::Component::Adaptor",

    has_planned_output => [
    ],
    has_provided_output => [
        fpkm_file => {
            is => 'Path',
        },
        tumor_sample_name => {
            is => 'Text',
        },
    ],
};

sub name {
    'fpkm';
}
