package Genome::Model::PhenotypeCorrelation::Command; 

use strict;
use warnings;

use Genome;

class Genome::Model::PhenotypeCorrelation::Command {
    is => 'Command::Tree',
    doc => 'operate on phenotype correlation models',
};

sub sub_command_category { 'type specific' }

1;
