package Genome::Model::SimpleAlignment;

use strict;
use warnings;

use Genome;

class Genome::Model::SimpleAlignment {
    is  => 'Genome::ModelDeprecated',
    has_param => [
       reference_sequence_name => { is => 'Text' },
    ],
};

sub _execute_build {
    die 'This model type can no longer be built.';
}

1;
