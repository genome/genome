package Genome::VariantReporting::Plan::FilterPlan;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Plan::FilterPlan {
    is => 'Genome::VariantReporting::Plan::Base',
};

sub category {
    'filters';
}

1;
