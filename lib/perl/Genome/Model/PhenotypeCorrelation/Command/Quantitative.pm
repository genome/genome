package Genome::Model::PhenotypeCorrelation::Command::Quantitative;

use Genome;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::Quantitative {
    is => 'Command',
    doc => 'commands for processing quantitative traits',
};

sub sub_command_category { 'trait specific' }
1;

