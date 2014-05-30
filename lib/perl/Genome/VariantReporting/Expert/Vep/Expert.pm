package Genome::VariantReporting::Expert::Vep::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Expert::Vep::Expert {
    is => 'Genome::VariantReporting::Expert::Base',
};

sub name {
    'vep';
}

sub priority {
    return 1;
}


1;
