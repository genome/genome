package Genome::Model::Tools::EpitopePrediction::Process;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::EpitopePrediction::Process {
    is => 'Genome::Process',
    has_input => [
        model => {
            is => 'Genome::Model::SomaticVariation',
        },
    ],
};

1;
