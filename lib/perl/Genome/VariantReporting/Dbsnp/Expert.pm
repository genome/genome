package Genome::VariantReporting::Dbsnp::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Dbsnp::Expert {
    is => 'Genome::VariantReporting::Framework::Component::Expert',
};

sub name {
    'dbsnp';
}


1;
