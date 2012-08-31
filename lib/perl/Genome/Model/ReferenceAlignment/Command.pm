use Genome;
use strict;
use warnings;

package Genome::Model::ReferenceAlignment::Command;

class Genome::Model::ReferenceAlignment::Command {
    is => 'Command',
    doc => 'operate on reference alignment models',
};

sub sub_command_category { 'type specific' }

1;

