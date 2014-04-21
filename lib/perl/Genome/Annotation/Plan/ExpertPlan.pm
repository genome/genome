package Genome::Annotation::Plan::ExpertPlan;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Plan::ExpertPlan {
    is => 'Genome::Annotation::Plan::Base',
};

sub category {
    'experts';
}

1;
