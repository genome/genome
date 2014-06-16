package Genome::VariantReporting::Expert::Dbsnp::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Expert::Dbsnp::Expert {
    is => 'Genome::VariantReporting::Expert::Base',
};

sub name {
    'dbsnp';
}


1;
