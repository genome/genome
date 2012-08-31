package Genome::Model::PhenotypeCorrelation::Command::Mendelian;

use Genome;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::Mendelian {
    is => 'Command',
    doc => 'commands for processing mendelian traits',
};

sub sub_command_category { 'trait specific' }
1;

