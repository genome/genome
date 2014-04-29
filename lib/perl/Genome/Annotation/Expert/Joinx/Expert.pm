package Genome::Annotation::Expert::Joinx::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::WorkflowBuilder::DAG;
use Genome::WorkflowBuilder::Command;

class Genome::Annotation::Expert::Joinx::Expert {
    is => 'Genome::Annotation::Expert::Base',
};

sub name {
    'joinx';
}


1;
