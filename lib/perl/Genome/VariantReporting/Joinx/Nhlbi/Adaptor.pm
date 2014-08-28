package Genome::VariantReporting::Joinx::Nhlbi::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Joinx::Nhlbi::Adaptor {
    is => 'Genome::VariantReporting::Joinx::Adaptor',
};

sub name {
    "nhlbi";
}

1;
