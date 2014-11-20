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
            doc => 'The sample to analyze',
        },
        fpkm_file => {
            is => 'Path',
            is_translated => 1,
            doc => 'The fpkm file to use for annotation',
        },
    ],
    doc => 'Annotate variants with FPKM value for the sample specified',
};

sub name {
    'fpkm';
}
