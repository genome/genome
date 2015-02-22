package Genome::Model::Tools::SvSim;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SvSim {
    is  => 'Command',
    is_abstract => 1,
};


sub help_brief {
    "Tools for sv simulation";
}

1;

