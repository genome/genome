package Genome::VariantReporting::Suite::Joinx::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Joinx::Expert {
    is => 'Genome::VariantReporting::Framework::Component::Expert',
    is_abstract => 1,
};

1;
