package Genome::VariantReporting::Suite::Joinx::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Joinx::Adaptor {
    is => 'Genome::VariantReporting::Framework::Component::Adaptor',
    is_abstract => 1,
    has_planned_output => [
        joinx_version => { is  => 'Version', },
        info_string => { is => 'Text', },
        vcf => {
            is => 'Path',
            is_translated => 1,
        }
    ],
};

1;
