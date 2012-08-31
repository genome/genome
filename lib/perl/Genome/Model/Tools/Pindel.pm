package Genome::Model::Tools::Pindel;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::Pindel {
    is => 'Command',
    has => [
    ],
};


sub help_brief {
    "Tools to run pindel or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
 gmt pindel ...    
EOS
}

sub help_detail {                           
    return <<EOS 
pindel stuff
EOS
}



1;

