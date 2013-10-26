package Genome::Model::PhenotypeCorrelation::Command::VcfAnnotation;

use strict;
use warnings;

use Genome;

class Genome::Model::PhenotypeCorrelation::Command::VcfAnnotation {
    is => 'Command::Tree',
};

sub sub_command_category { 'type specific' }

1;
