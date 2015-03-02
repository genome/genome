package Genome::VariantReporting::Suite::Joinx::Dbsnp::Run;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Suite::Joinx::Dbsnp::Run {
    is => 'Genome::VariantReporting::Suite::Joinx::Run',
    doc => 'Annotate vcf with results from dbsnp vcf',
};

sub name {
    return 'dbsnp';
}

1;
