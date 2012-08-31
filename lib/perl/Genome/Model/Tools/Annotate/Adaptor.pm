package Genome::Model::Tools::Annotate::Adaptor;

use strict;
use warnings;

use Genome;     

class Genome::Model::Tools::Annotate::Adaptor {
    is => 'Command',
};

sub sub_command_sort_position { 15 }

sub help_brief {
    "tools used for adapting various file formats into a format the annotation tool can accept"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools annotate adaptor ...    
EOS
}

