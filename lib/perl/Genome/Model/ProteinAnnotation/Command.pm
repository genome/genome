package Genome::Model::ProteinAnnotation::Command;

use strict;
use warnings;

use Genome;

class Genome::Model::ProteinAnnotation::Command {
    is => 'Command::Tree',
    doc => 'operate on protein annotation models',
};

sub sub_command_category { 'type specific' }

1;
