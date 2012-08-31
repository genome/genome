package Genome::Model::Tools::Assembly;

use strict;
use warnings;

use Genome;            

class Genome::Model::Tools::Assembly {
    is => 'Command',
    has => [ ],
};

#sub sub_command_sort_position { 15 }

sub help_brief {
    'Tools for working with pcap assemblies.'
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools assembly ...
EOS
}

sub xhelp_detail {                           
    return <<EOS 
EOS
}

1;
