package Genome::VariantReporting::Expert::Joinx::Dbsnp::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Expert::Joinx::Dbsnp::Expert {
    is => 'Genome::VariantReporting::Expert::Joinx::Expert',
};

sub name {
    'dbsnp';
}

1;
