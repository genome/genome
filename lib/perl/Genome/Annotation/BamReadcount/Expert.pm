package Genome::Annotation::BamReadcount::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::WorkflowBuilder::DAG;
use Genome::WorkflowBuilder::Command;

class Genome::Annotation::BamReadcount::Expert {
    is => 'Genome::Annotation::ExpertBase',
};

sub name {
    'bam-readcount';
}

sub dag {
    return;
}

1;
