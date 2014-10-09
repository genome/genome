package Genome::VariantReporting::Suite::Joinx::Dbsnp::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Joinx::Dbsnp::Expert {
    is => 'Genome::VariantReporting::Suite::Joinx::Expert',
};

sub name {
    'dbsnp';
}

1;
