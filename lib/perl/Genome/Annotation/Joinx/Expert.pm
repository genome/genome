package Genome::Annotation::Joinx::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::WorkflowBuilder::DAG;
use Genome::WorkflowBuilder::Command;

class Genome::Annotation::Joinx::Expert {
    is => 'Genome::Annotation::ExpertBase',
};

sub name {
    'joinx';
}


1;
