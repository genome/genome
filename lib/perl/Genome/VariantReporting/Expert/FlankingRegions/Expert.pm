package Genome::VariantReporting::Expert::FlankingRegions::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Expert::FlankingRegions::Expert {
    is => 'Genome::VariantReporting::Framework::Component::Expert',
};

sub name {
    'flanking-regions';
}

1;
