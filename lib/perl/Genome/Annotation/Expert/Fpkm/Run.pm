package Genome::Annotation::Expert::Fpkm::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Expert::Fpkm::Run {
    is => 'Genome::Annotation::Expert::CommandBase',
    has_input => [
        fpkm_file => {is => 'Path'},
        tumor_sample_name => {is => 'Text'},
    ],
};

sub name {
    'fpkm';
}

sub result_class {
    'Genome::Annotation::Expert::Fpkm::RunResult';
}
