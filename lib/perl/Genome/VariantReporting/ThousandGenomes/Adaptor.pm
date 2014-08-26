package Genome::VariantReporting::ThousandGenomes::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::ThousandGenomes::Adaptor {
    is => 'Genome::VariantReporting::Joinx::Adaptor',
};

sub name {
    "1kg";
}

1;
