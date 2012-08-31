package Genome::Model::Tools::Library;

use strict;
use warnings;

use Genome;            

class Genome::Model::Tools::Library {
    is => 'Command',
    has => [ ],
};

sub sub_command_sort_position { 15 }

sub help_brief {
    "Tools for working with a model's library."
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools library ...
EOS
}

sub xhelp_detail {                           
    return <<EOS 
EOS
}

1;

