package Genome::VariantReporting::Joinx::HomoPolymer::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Joinx::HomoPolymer::Expert {
    is => 'Genome::VariantReporting::Framework::Component::Expert',
};

sub name {
    'homo-polymer';
}

1;
