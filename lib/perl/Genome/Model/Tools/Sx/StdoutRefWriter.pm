package Genome::Model::Tools::Sx::StdoutRefWriter;

use strict;
use warnings;

require Storable;

class Genome::Model::Tools::Sx::StdoutRefWriter { 
};

sub write {
    Storable::store_fd($_[1], \*STDOUT);
    return 1;
}

sub flush { return 1; }

1;

