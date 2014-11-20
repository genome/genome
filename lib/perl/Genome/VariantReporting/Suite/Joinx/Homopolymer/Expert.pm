package Genome::VariantReporting::Suite::Joinx::Homopolymer::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Joinx::Homopolymer::Expert {
    is => 'Genome::VariantReporting::Framework::Component::Expert',
};

sub name {
    'homopolymer';
}

1;
