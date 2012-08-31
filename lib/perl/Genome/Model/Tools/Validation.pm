package Genome::Model::Tools::Validation;

use strict;
use warnings;

use Genome;            

class Genome::Model::Tools::Validation {
    is => 'Command',
    has => [ ],
};

sub sub_command_sort_position { 15 }

sub help_brief {
    "Tools for working with validation results"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt validation ...
EOS
}

sub xhelp_detail {                           
    return <<EOS 
EOS
}

1;

