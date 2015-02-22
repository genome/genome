package Genome::VariantReporting::Suite::Joinx::Dbsnp::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Joinx::Dbsnp::Adaptor {
    is => 'Genome::VariantReporting::Suite::Joinx::Adaptor',
    doc => 'Annotate vcf with results from dbsnp vcf',
};

sub name {
    "dbsnp";
}

1;
