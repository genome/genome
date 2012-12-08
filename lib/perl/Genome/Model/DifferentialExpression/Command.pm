package Genome::Model::DifferentialExpression::Command;

use strict;
use warnings;

use Genome;

class Genome::Model::DifferentialExpression::Command {
    is => 'Command::Tree',
    doc => 'operate on differential expression models',
};

sub sub_command_category { 'type specific' }

1;

