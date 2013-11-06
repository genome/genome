package Genome::Model::RnaSeq::Command::DetectFusions::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions::Base {
    is => 'Genome::Command::Base',
    is_abstract => 1,
    has_input => [
        detector_version => {
            is => 'Text',
            doc => 'the version of the fusion detector to run',
        },
        detector_params => {
            is => 'Text',
            doc => 'parameters for the chosen fusion detector',
        },
        build => {
            is => "Genome::Model::Build::RnaSeq",
        },
    ],
    has_optional_output => [
        result => {
            is => 'Genome::SoftwareResult',
        },
    ],
};
