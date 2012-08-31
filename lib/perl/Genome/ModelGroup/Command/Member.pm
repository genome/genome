package Genome::ModelGroup::Command::Member;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::Command::Member {
    is => 'Command::Tree',
    doc => "work with the members of model-groups",
};

sub help_synopsis {
    return <<"EOS"
    work with the members of model-groups   
EOS
}

sub help_brief {
    return "work with the members of model-groups";
}

1;

