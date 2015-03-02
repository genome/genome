package Genome::VariantReporting::Suite::Joinx::ThousandGenomes::Run;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Suite::Joinx::ThousandGenomes::Run {
    is => 'Genome::VariantReporting::Suite::Joinx::Run',
    doc => 'Annotate vcf with results from 1000Genomes vcf',
};

sub name {
    return '1kg';
}

1;
