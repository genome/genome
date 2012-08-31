package Genome::Model::Tools::Graph;

use strict;
use warnings;

use Genome;     

class Genome::Model::Tools::Graph {
    is => 'Command',
};

sub sub_command_sort_position { 14 }

sub help_brief {
    "graph/image generators"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools graph ...    
EOS
}

1;

