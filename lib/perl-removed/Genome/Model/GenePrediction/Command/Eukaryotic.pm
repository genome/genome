package Genome::Model::GenePrediction::Command::Eukaryotic;

use Genome;
use Carp 'confess';

use strict;
use warnings;

class Genome::Model::GenePrediction::Command::Eukaryotic {
    is => 'Command',
    doc => 'operate on gene prediction models',
};

sub sub_command_category { 'type specific' }

1;
