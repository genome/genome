package Genome::VariantReporting::Dbsnp::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Dbsnp::Adaptor {
    is => 'Genome::VariantReporting::Joinx::Adaptor',
};

sub name {
    "dbsnp";
}

1;
