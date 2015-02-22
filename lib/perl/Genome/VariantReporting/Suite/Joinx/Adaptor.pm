package Genome::VariantReporting::Suite::Joinx::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Joinx::Adaptor {
    is => 'Genome::VariantReporting::Framework::Component::Adaptor',
    is_abstract => 1,
    has_planned_output => [
        joinx_version => {
            is  => 'Version',
            doc => "joinx version to be used.",
        },
        info_string => {
            is => 'Text',
            doc => 'Field ids to embed from the annotation VCF. Use colons to separate multiple field descriptors.',
        },
        vcf => {
            is => 'Path',
            is_translated => 1,
            doc => 'Vcf File containing annotation',
        }
    ],
};

1;
