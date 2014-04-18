package Genome::Annotation::Vep::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::WorkflowBuilder::DAG;
use Genome::WorkflowBuilder::Command;

class Genome::Annotation::Vep::Expert {
    is => 'Genome::Annotation::ExpertBase',
};

sub name {
    'vep';
}

sub dag {
    return;
}

1;
