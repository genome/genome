package Genome::Annotation::Plan::InterpreterPlan;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Plan::InterpreterPlan {
    is => 'Genome::Annotation::Plan::Base',
};

sub category {
    'interpreters';
}

1;
