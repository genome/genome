package Genome::VariantReporting::Dbsnp::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Dbsnp::Adaptor {
    is => 'Genome::VariantReporting::Framework::Component::Adaptor',

    has_planned_output => [
        joinx_version => { is  => 'Version', },
        info_string => { is => 'Text', },
    ],
    has_provided_output => [
        dbsnp_vcf => {
            is => 'Path',
        }
    ],
};

sub name {
    "dbsnp";
}

1;
