package Genome::VariantReporting::Framework::Plan::FilterPlan;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Framework::Plan::FilterPlan {
    is => 'Genome::VariantReporting::Framework::Plan::Base',
};

sub category {
    'filters';
}

1;
