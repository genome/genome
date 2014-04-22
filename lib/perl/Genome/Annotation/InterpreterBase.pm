package Genome::Annotation::InterpreterBase;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::InterpreterBase {
    is => 'Genome::Annotation::ComponentBase',
    is_abstract => 1,
};

sub requires_experts {
    return ();
}

1;
