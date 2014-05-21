package Genome::Annotation::Expert::Dbsnp::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Expert::Dbsnp::Expert {
    is => 'Genome::Annotation::Expert::Base',
};

sub name {
    'dbsnp';
}


1;
