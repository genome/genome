package Genome::VariantReporting::Joinx::Dbsnp::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Joinx::Dbsnp::Adaptor {
    is => 'Genome::VariantReporting::Joinx::Adaptor',
};

sub name {
    "dbsnp";
}

1;
