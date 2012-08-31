package Genome::Model::Tools::Array;

use strict;
use warnings;

use Genome;            

class Genome::Model::Tools::Array {
    is => 'Command',
    has => [ ],
};

sub sub_command_sort_position { 15 }

sub help_brief {
    'Tools for working with microarray data of various kinds.'
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools array ...
EOS
}

sub xhelp_detail {                           
    return <<EOS 
EOS
}

1;

