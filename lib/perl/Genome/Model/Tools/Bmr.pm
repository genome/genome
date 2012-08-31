package Genome::Model::Tools::Bmr;

use strict;
use warnings;

use Genome;            

class Genome::Model::Tools::Bmr {
    is => 'Command',
    has => [ ],
};

sub sub_command_sort_position { 15 }

sub help_brief {
    'Tools for working with background mutation rates.'
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt bmr ...
EOS
}

sub xhelp_detail {                           
    return <<EOS 
EOS
}

1;

