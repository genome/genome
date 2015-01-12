package Genome::Model::PhenotypeCorrelation::Command::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::PhenotypeCorrelation::Command::Base {
    is => "Command::V2",
    is_abstract => 1,
    has_input => [
        build => {
            is => 'Genome::Model::Build::PhenotypeCorrelation',
        },
    ],
    has => [
        processing_profile => {
            is => 'Genome::ProcessingProfile::PhenotypeCorrelation',
            via => 'build',
            to => 'processing_profile',
        },
    ],
    doc => 'base class for the delegate classes for PhenotypeCorrelation builds',
};

1;
