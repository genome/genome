package Genome::Annotation::Plan::FilterPlan;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Plan::FilterPlan {
    is => 'Genome::Annotation::Plan::Base',
};

sub category {
    'filters';
}

1;
