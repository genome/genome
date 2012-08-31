package Genome::Model::PhenotypeCorrelation::Command; 

use strict;
use warnings;

use Genome;

class Genome::Model::PhenotypeCorrelation::Command {
    is => 'Command::Tree',
};

sub sub_command_category { 'type specific' }

1;
