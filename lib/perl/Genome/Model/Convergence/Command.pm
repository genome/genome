use Genome;
use strict;
use warnings;

package Genome::Model::Convergence::Command;

class Genome::Model::Convergence::Command {
    is => 'Command',
    doc => 'operate on convergence models',
};

sub sub_command_category { 'type specific' }

1;

