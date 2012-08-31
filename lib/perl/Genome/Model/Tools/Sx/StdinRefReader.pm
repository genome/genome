package Genome::Model::Tools::Sx::StdinRefReader;

use strict;
use warnings;

require Storable;

class Genome::Model::Tools::Sx::StdinRefReader { 
}; 

sub read {
    my $self = shift;

    my $ref = eval { Storable::fd_retrieve(*STDIN) };

    return $ref;
}

1;

