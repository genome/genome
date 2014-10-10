package Genome::VariantReporting::Suite::Fpkm::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;
use File::Spec;

class Genome::VariantReporting::Suite::Fpkm::Adaptor {
    is => "Genome::VariantReporting::Framework::Component::Adaptor",
    has_planned_output => [
        sample_name => {
            is => 'Text',
            is_translated => 1,
        }
    ],
    has_provided_output => [
        fpkm_file => {
            is => 'Path',
        },
    ],
};

sub name {
    'fpkm';
}
