package Genome::VariantReporting::Plan::InterpreterPlan;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Plan::InterpreterPlan {
    is => 'Genome::VariantReporting::Plan::Base',
};

sub category {
    'interpreters';
}

1;
