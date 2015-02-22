package Genome::VariantReporting::Suite::Joinx::ThousandGenomes::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Joinx::ThousandGenomes::Adaptor {
    is => 'Genome::VariantReporting::Suite::Joinx::Adaptor',
    doc => 'Annotate vcf with results from 1000Genomes vcf',
};

sub name {
    "1kg";
}

1;
