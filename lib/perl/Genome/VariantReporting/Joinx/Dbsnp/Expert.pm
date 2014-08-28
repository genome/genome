package Genome::VariantReporting::Joinx::Dbsnp::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Joinx::Dbsnp::Expert {
    is => 'Genome::VariantReporting::Joinx::Expert',
};

sub name {
    'dbsnp';
}

1;
