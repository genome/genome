package Genome::VariantReporting::Vep::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Vep::Expert {
    is => 'Genome::VariantReporting::Component::Expert',
};

sub name {
    'vep';
}

sub priority {
    return 1;
}


1;
