package Genome::Model::Tools::Samtools;

use strict;
use warnings;

use Genome;                         # >above< ensures YOUR copy is used during development

class Genome::Model::Tools::Samtools {
    is => 'Command',
};

sub sub_command_sort_position { 14 }

sub help_brief {
    "Tools to run Samtools or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools Samtools ...    
EOS
}

sub help_detail {                           
    return <<EOS 
EOS
}

1;

