package Genome::Model::Tools::EpitopePrediction::Process;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::EpitopePrediction::Process {
    is => 'Genome::Process',
    has_input => [
        build => {
            is => 'Genome::Model::Build::SomaticVariation',
        },
    ],
};

1;
