use Genome;
use strict;
use warnings;

package Genome::Model::ReferenceSequence::Command;

class Genome::Model::ReferenceSequence::Command {
    is => 'Command',
    doc => 'operate on reference alignment models',
};

sub sub_command_category { 'type specific' }

1;

