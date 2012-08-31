use Genome;
use strict;
use warnings;

package Genome::Model::Somatic::Command;

class Genome::Model::Somatic::Command {
    is => 'Command',
    doc => 'operate on somatic models',
};

sub sub_command_category { 'type specific' }

1;

