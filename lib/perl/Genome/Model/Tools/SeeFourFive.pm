package Genome::Model::Tools::SeeFourFive;

use strict;
use warnings;

use Genome;            

class Genome::Model::Tools::SeeFourFive {
    is => 'Command',
    has => [ ],
};

sub sub_command_sort_position { 16 }

sub help_brief {
    'Tools for working with C4.5 and SeeFive.'
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools c5
EOS
}

1;

