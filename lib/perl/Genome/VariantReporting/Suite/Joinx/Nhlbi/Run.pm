package Genome::VariantReporting::Suite::Joinx::Nhlbi::Run;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Suite::Joinx::Nhlbi::Run {
    is => 'Genome::VariantReporting::Suite::Joinx::Run',
    doc => 'Annotate vcf with results from NHLBI vcf',
};

sub name {
    return 'nhlbi';
}

1;
