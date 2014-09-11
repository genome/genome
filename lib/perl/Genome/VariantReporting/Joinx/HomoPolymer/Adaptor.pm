package Genome::VariantReporting::Joinx::HomoPolymer::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Joinx::HomoPolymer::Adaptor {
    is => "Genome::VariantReporting::Framework::Component::Adaptor",

    has_planned_output => [
        joinx_version => {
            is => 'Version',
        },
        max_length => {
            is => 'Integer',
        },
        info_string => { 
            is => 'Text', 
        },
    ],
    has_provided_output => [
        homopolymer_list_id => {
            is => 'String',
        },
    ],
};

sub name {
    'homo-polymer';
}

1;
