package Genome::Model::Tools::Consed;

use strict;
use warnings;

use Genome;                         # >above< ensures YOUR copy is used during development

class Genome::Model::Tools::Consed {
    is => 'Command',
};

sub sub_command_sort_position { 13 }

sub help_brief {
    "Tools to run consed or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools consed ...    
EOS
}

sub help_detail {                           
    return <<EOS 
EOS
}

1;

