package Genome::Model::GenePrediction::Command;

use Genome;
use Carp 'confess';

use strict;
use warnings;

class Genome::Model::GenePrediction::Command {
    is => 'Command',
    doc => 'operate on gene prediction models',
};

sub sub_command_category { 'type specific' }

1;
