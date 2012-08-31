package Genome::Model::Tools::FastTier;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::FastTier {
    is => 'Command',
    has => [
    ],
};


sub help_brief {
    "Tools to run maq or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
 gmt fast-tier ...    
EOS
}

sub help_detail {                           
    return <<EOS 
fast tiering stuff
EOS
}



1;

