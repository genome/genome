package Genome::Model::Tools::Assembly::Repair;

use strict;
use warnings;

use Genome;            

class Genome::Model::Tools::Assembly::Repair {
    is => 'Command',
    has => [ ],
};

#sub sub_command_sort_position { 15 }

sub help_brief {
    'Tools for repair and improvement of assembly data.'
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools assembly repair ...
EOS
}

sub xhelp_detail {                           
    return <<EOS 
EOS
}

1;

