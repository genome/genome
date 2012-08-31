package Genome::Model::Tools::Import;

use strict;
use warnings;

use Genome;            

class Genome::Model::Tools::Import {
    is => 'Command',
    has => [ ],
};

sub sub_command_sort_position { 15 }

sub help_brief {
    'Tools for working with SNP files of various kinds.'
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools import ...
EOS
}

sub xhelp_detail {                           
    return <<EOS 
EOS
}

1;

