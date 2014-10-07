package Genome::VariantReporting::Expert::Joinx::Dbsnp::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Expert::Joinx::Dbsnp::Adaptor {
    is => 'Genome::VariantReporting::Expert::Joinx::Adaptor',
};

sub name {
    "dbsnp";
}

1;
