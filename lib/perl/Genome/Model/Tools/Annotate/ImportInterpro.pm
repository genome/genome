package Genome::Model::Tools::Annotate::ImportInterpro;

use strict;
use warnings;

use Genome;     

class Genome::Model::Tools::Annotate::ImportInterpro {
    is => 'Command',
};

sub sub_command_sort_position { 15 }

sub help_brief {
    "tools used for importing interpro results"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools annotate import-interpro ...    
EOS
}

