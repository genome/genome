package Genome::Model::Tools::Soap;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Soap {
    is => 'Command::Tree',
};

sub help_detail {
    return 'De novo assemble and align with soap';
}

1;
