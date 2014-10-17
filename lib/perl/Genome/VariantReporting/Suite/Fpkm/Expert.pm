package Genome::VariantReporting::Suite::Fpkm::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Fpkm::Expert {
    is => 'Genome::VariantReporting::Framework::Component::Expert',
};

sub name {
    'fpkm';
}


1;
