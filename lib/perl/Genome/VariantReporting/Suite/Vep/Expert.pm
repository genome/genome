package Genome::VariantReporting::Suite::Vep::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Vep::Expert {
    is => 'Genome::VariantReporting::Framework::Component::Expert',
};

sub name {
    'vep';
}

sub priority {
    return 1;
}


1;
