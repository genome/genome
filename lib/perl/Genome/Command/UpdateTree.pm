package Genome::Command::UpdateTree;

use strict;
use warnings;

use Genome;
      
class Genome::Command::UpdateTree {
    is => 'Command::Tree',
    doc => 'CRUD update tree class.',
};

sub sub_command_sort_position { .3 };

1;

