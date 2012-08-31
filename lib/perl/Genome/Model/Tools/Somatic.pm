package Genome::Model::Tools::Somatic;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::Somatic {
    is => 'Command',
    has => [
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to run maq or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
 gmt somatic ...    
EOS
}

sub help_detail {                           
    return <<EOS 
More information about the maq suite of tools can be found at http://maq.sourceforege.net.
EOS
}



1;

