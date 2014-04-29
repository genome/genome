package Genome::Annotation::Interpreter::Base;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Interpreter::Base {
    is => 'Genome::Annotation::ComponentBase',
    is_abstract => 1,
};

sub requires_experts {
    return ();
}

1;
