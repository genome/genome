package Genome::ModelGroup::Command;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::Command {
    is => 'Command::Tree',
    doc => "work with model-groups",
};

sub help_synopsis {
    return 'genome model-group' 
}

sub help_detail {                           
    return 'Top level command to hold commands for working with model-groups'
}

1;
