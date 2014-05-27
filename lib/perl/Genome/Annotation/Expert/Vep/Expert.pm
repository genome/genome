package Genome::Annotation::Expert::Vep::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Expert::Vep::Expert {
    is => 'Genome::Annotation::Expert::Base',
};

sub name {
    'vep';
}

sub priority {
    return 1;
}


1;
