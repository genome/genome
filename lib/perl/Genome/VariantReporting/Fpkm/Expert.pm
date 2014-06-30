package Genome::VariantReporting::Fpkm::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Fpkm::Expert {
    is => 'Genome::VariantReporting::Framework::Component::Expert',
};

sub name {
    'fpkm';
}


1;
