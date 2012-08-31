use Genome;
use strict;
use warnings;

package Genome::Model::SomaticValidation::Command;

class Genome::Model::SomaticValidation::Command {
    is => 'Command',
    doc => 'operate on somatic-validation models',
};

sub sub_command_category { 'type specific' }

1;

