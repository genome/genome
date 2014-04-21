package Genome::Annotation::FilterBase;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::FilterBase {
    is => 'Genome::Annotation::ComponentBase',
    is_abstract => 1,
};

1;
