package Genome::Model::PhenotypeCorrelation::Command::CaseControl;

use Genome;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::CaseControl {
    is => 'Command',
    doc => 'commands for processing case/control type traits',
};

sub sub_command_category { 'trait specific' }
1;

