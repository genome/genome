package Genome::Annotation::ReporterBase;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::ReporterBase {
    is => 'Genome::Annotation::ComponentBase',
    is_abstract => 1,
};

sub name {
    die "abstract";
}

1;
