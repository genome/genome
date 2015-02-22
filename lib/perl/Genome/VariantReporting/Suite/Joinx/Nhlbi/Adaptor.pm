package Genome::VariantReporting::Suite::Joinx::Nhlbi::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Joinx::Nhlbi::Adaptor {
    is => 'Genome::VariantReporting::Suite::Joinx::Adaptor',
    doc => 'Annotate vcf with results from NHLBI vcf',
};

sub name {
    "nhlbi";
}

1;
