package Genome::Model::Tools::Velvet;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Velvet {
    is => 'Command::Tree',
};

sub help_detail {
    return 'De novo assemble with velvet';
}

1;
