package Genome::VariantReporting::Nhlbi::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Nhlbi::Adaptor {
    is => 'Genome::VariantReporting::Joinx::Adaptor',
};

sub name {
    "nhlbi";
}

1;
