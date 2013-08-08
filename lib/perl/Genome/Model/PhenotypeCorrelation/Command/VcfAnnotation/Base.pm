package Genome::Model::PhenotypeCorrelation::Command::VcfAnnotation::Base;

use Carp qw/confess/;
use Genome;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::VcfAnnotation::Base {
    is => "Genome::Command::Base",
    has_input => [
        input_file => {
            is => "Text",
            doc => "Input vcf file",
        },
        output_file => {
            is => "Text",
            is_output => 1,
            doc => "Output vcf file",
        },
        species_name => {
            is => "Text",
            default_value => "human",
        }
    ]
};

1;
