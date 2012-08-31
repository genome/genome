use Genome;
use strict;
use warnings;

package Genome::Model::SomaticVariation::Command;

class Genome::Model::SomaticVariation::Command {
    is => 'Command',
    doc => 'operate on somatic-variation models',
};

sub sub_command_category { 'type specific' }

1;

