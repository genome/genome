package Genome::VariantReporting::Nhlbi::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Nhlbi::Expert {
    is => 'Genome::VariantReporting::Joinx::Expert',
};

sub name {
    'nhlbi';
}

1;
