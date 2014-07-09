package Genome::VariantReporting::Framework::Test::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Framework::Test::Expert {
    is => 'Genome::VariantReporting::Framework::Component::Expert',
};

sub name {
    '__test__';
}


1;
