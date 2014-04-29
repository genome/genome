package Genome::Annotation::Expert::Vep::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::WorkflowBuilder::DAG;
use Genome::WorkflowBuilder::Command;

class Genome::Annotation::Expert::Vep::Expert {
    is => 'Genome::Annotation::Expert::Base',
};

sub name {
    'vep';
}

sub dag {
    return;
}

1;
