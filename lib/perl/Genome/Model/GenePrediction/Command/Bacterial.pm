package Genome::Model::GenePrediction::Command::Bacterial;

use Genome;

use strict;
use warnings;

class Genome::Model::GenePrediction::Command::Bacterial {
    is => 'Command',
    doc => 'operate on gene prediction models',
};

sub sub_command_category { 'type specific' }
1;

