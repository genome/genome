package Genome::Annotation::Expert::Fpkm::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Expert::Fpkm::Expert {
    is => 'Genome::Annotation::Expert::Base',
};

sub name {
    'fpkm';
}


1;
