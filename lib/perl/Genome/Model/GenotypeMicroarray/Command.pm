use Genome;
use strict;
use warnings;

package Genome::Model::GenotypeMicroarray::Command;

class Genome::Model::GenotypeMicroarray::Command {
    is => 'Command',
    doc => 'operate on genotype microarray models',
};

sub sub_command_category { 'type specific' }

1;

