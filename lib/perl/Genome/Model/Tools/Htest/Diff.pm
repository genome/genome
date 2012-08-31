
package Genome::Model::Tools::Htest::Diff;

use strict;
use warnings;

use Genome;
use Command; 

class Genome::Model::Tools::Htest::Diff {
    is => 'Command',
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "hypothesize and test one or more sequence variations"
}

1;

