use Genome;

use strict;
use warnings;

package Genome::Model::SomaticCapture::Command;

class Genome::Model::SomaticCapture::Command {
    is => 'Command',
    doc => 'operate on somatic-capture models',
};

sub sub_command_category { 'type specific' }

1;
